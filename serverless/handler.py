"""RunPod serverless handler for DockQ.

Wraps the installed `DockQ` command line tool; DockQ itself is unchanged.

Job input (``job["input"]``), for each of ``model`` and ``native`` give ONE of:
    "<name>":        structure file contents as plain text (PDB or mmCIF)
    "<name>_base64": base64-encoded file (plain or gzip-compressed)
    "<name>_url":    http(s) URL to download the file from (plain or .gz)

Optional DockQ options (same meaning as the CLI flags):
    "mapping": "AB:HL", "small_molecule": false, "capri_peptide": false,
    "no_align": false, "allowed_mismatches": 0, "n_cpu": <cpus>, "max_chunk": 512

Returns the same dictionary DockQ writes with ``--json``.
"""

import base64
import json
import os
import subprocess
import tempfile
import urllib.request

import runpod

GZIP_MAGIC = b"\x1f\x8b"
TIMEOUT = int(os.environ.get("DOCKQ_TIMEOUT", "3600"))
MAX_DOWNLOAD_BYTES = int(os.environ.get("DOCKQ_MAX_DOWNLOAD_BYTES", str(500 * 1024**2)))

BOOL_FLAGS = ("small_molecule", "capri_peptide", "no_align")
INT_OPTIONS = ("allowed_mismatches", "n_cpu", "max_chunk")


def _get_bytes(job_input, name):
    given = [k for k in (name, f"{name}_base64", f"{name}_url") if job_input.get(k)]
    if len(given) != 1:
        raise ValueError(
            f"Provide exactly one of '{name}', '{name}_base64' or '{name}_url'"
        )
    key = given[0]
    value = job_input[key]
    if key == name:
        return value.encode()
    if key.endswith("_base64"):
        return base64.b64decode(value)
    if not value.startswith(("http://", "https://")):
        raise ValueError(f"'{key}' must be an http(s) URL")
    with urllib.request.urlopen(value, timeout=300) as resp:
        data = resp.read(MAX_DOWNLOAD_BYTES + 1)
    if len(data) > MAX_DOWNLOAD_BYTES:
        raise ValueError(f"'{key}' is larger than {MAX_DOWNLOAD_BYTES} bytes")
    return data


def _write_structure(workdir, name, data):
    # DockQ picks gzip.open based on the ".gz" suffix and auto-detects PDB/mmCIF
    suffix = ".gz" if data[:2] == GZIP_MAGIC else ""
    path = os.path.join(workdir, f"{name}{suffix}")
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def _build_command(job_input, model_path, native_path, json_path):
    cmd = ["DockQ", model_path, native_path, "--json", json_path, "--short"]
    for flag in BOOL_FLAGS:
        if job_input.get(flag):
            cmd.append(f"--{flag}")
    for opt in INT_OPTIONS:
        if job_input.get(opt) is not None:
            cmd += [f"--{opt}", str(int(job_input[opt]))]
    if "n_cpu" not in job_input:
        cmd += ["--n_cpu", str(os.cpu_count() or 1)]
    if job_input.get("mapping"):
        cmd += ["--mapping", str(job_input["mapping"])]
    return cmd


def handler(job):
    job_input = job.get("input") or {}
    try:
        with tempfile.TemporaryDirectory() as workdir:
            model_path = _write_structure(workdir, "model", _get_bytes(job_input, "model"))
            native_path = _write_structure(workdir, "native", _get_bytes(job_input, "native"))
            json_path = os.path.join(workdir, "result.json")

            cmd = _build_command(job_input, model_path, native_path, json_path)
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
            if proc.returncode != 0 or not os.path.exists(json_path):
                return {
                    "error": "DockQ failed",
                    "returncode": proc.returncode,
                    "stderr": proc.stderr[-4000:],
                    "stdout": proc.stdout[-4000:],
                }

            with open(json_path) as fh:
                result = json.load(fh)
    except (ValueError, OSError, subprocess.TimeoutExpired) as exc:
        return {"error": str(exc)}

    # Temp paths are meaningless to the caller
    result.pop("model", None)
    result.pop("native", None)
    return result


if __name__ == "__main__":
    runpod.serverless.start({"handler": handler})
