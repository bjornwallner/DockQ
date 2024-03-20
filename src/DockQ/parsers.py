import numpy as np
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBParser import as_handle
import Bio
import warnings

class PDBParser(Bio.PDB.PDBParser):
    def get_structure(self, id, file, chains):
        """Return the structure.

        Arguments:
         - id - string, the id that will be used for the structure
         - file - name of the PDB file OR an open filehandle

        """
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            self.header = None
            self.trailer = None
            # Make a StructureBuilder instance (pass id of structure as parameter)
            self.structure_builder.init_structure(id)

            with as_handle(file) as handle:
                lines = handle.readlines()
                if not lines:
                    raise ValueError("Empty file.")
                self._parse(lines, chains)

            self.structure_builder.set_header(self.header)
            # Return the Structure instance
            structure = self.structure_builder.get_structure()

        return structure


    def _parse(self, header_coords_trailer, chains):
        """Parse the PDB file (PRIVATE)."""
        # Extract the header; return the rest of the file
        self.header, coords_trailer = self._get_header(header_coords_trailer)
        # Parse the atomic data; return the PDB file trailer
        self.trailer = self._parse_coordinates(coords_trailer, chains)


    def _parse_coordinates(self, coords_trailer, chains=[]):
        """Parse the atomic data in the PDB file (PRIVATE)."""
        allowed_records = {
            "ATOM  ",
            "HETATM",
            "MODEL ",
            "ENDMDL",
            "TER   ",
        }
        local_line_counter = 0
        structure_builder = self.structure_builder
        current_model_id = 0
        # Flag we have an open model
        model_open = 0
        current_chain_id = None
        current_segid = None
        current_residue_id = None
        current_resname = None

        for i in range(len(coords_trailer)):
            line = coords_trailer[i].rstrip("\n")
            record_type = line[0:6]
            global_line_counter = self.line_counter + local_line_counter + 1
            structure_builder.set_line_counter(global_line_counter)
            if not line.strip():
                continue  # skip empty lines
            elif record_type == "HETATM":
                continue
            elif record_type == "ATOM  ":
                # Initialize the Model - there was no explicit MODEL record
                if not model_open:
                    structure_builder.init_model(current_model_id)
                    current_model_id += 1
                    model_open = 1
                chainid = line[21]
                if chains and chainid not in chains:
                    continue
                fullname = line[12:16]
                # get rid of whitespace in atom names
                split_list = fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    name = fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    name = split_list[0]
                altloc = line[16]
                resname = line[17:20].strip()
                try:
                    serial_number = int(line[6:11])
                except Exception:
                    serial_number = 0
                resseq = int(line[22:26].split()[0])  # sequence identifier
                icode = line[26]  # insertion code
                if record_type == "HETATM":  # hetero atom flag
                    if resname == "HOH" or resname == "WAT":
                        hetero_flag = "W"
                    else:
                        hetero_flag = "H"
                else:
                    hetero_flag = " "
                residue_id = (hetero_flag, resseq, icode)
                # atomic coordinates
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    # Should we allow parsing to continue in permissive mode?
                    # If so, what coordinates should we default to?  Easier to abort!
                    raise PDBConstructionException(
                        "Invalid or missing coordinate(s) at line %i."
                        % global_line_counter
                    ) from None
                coord = np.array((x, y, z), "f")

                # occupancy & B factor
                if not self.is_pqr:
                    try:
                        occupancy = float(line[54:60])
                    except Exception:
                        self._handle_PDB_exception(
                            "Invalid or missing occupancy", global_line_counter
                        )
                        occupancy = None  # Rather than arbitrary zero or one
                    if occupancy is not None and occupancy < 0:
                        # TODO - Should this be an error in strict mode?
                        # self._handle_PDB_exception("Negative occupancy",
                        #                            global_line_counter)
                        # This uses fixed text so the warning occurs once only:
                        warnings.warn(
                            "Negative occupancy in one or more atoms",
                            PDBConstructionWarning,
                        )
                    try:
                        bfactor = float(line[60:66])
                    except Exception:
                        self._handle_PDB_exception(
                            "Invalid or missing B factor", global_line_counter
                        )
                        bfactor = 0.0  # PDB uses a default of zero if missing

                elif self.is_pqr:
                    # Attempt to parse charge and radius fields
                    try:
                        pqr_charge = float(line[54:62])
                    except Exception:
                        self._handle_PDB_exception(
                            "Invalid or missing charge", global_line_counter
                        )
                        pqr_charge = None  # Rather than arbitrary zero or one
                    try:
                        radius = float(line[62:70])
                    except Exception:
                        self._handle_PDB_exception(
                            "Invalid or missing radius", global_line_counter
                        )
                        radius = None
                    if radius is not None and radius < 0:
                        # In permissive mode raise fatal exception.
                        message = "Negative atom radius"
                        self._handle_PDB_exception(message, global_line_counter)
                        radius = None

                segid = line[72:76]
                element = line[76:78].strip().upper()
                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)

                if not self.is_pqr:
                    # init atom with pdb fields
                    try:
                        structure_builder.init_atom(
                            name,
                            coord,
                            bfactor,
                            occupancy,
                            altloc,
                            fullname,
                            serial_number,
                            element,
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif self.is_pqr:
                    try:
                        structure_builder.init_atom(
                            name,
                            coord,
                            pqr_charge,
                            radius,
                            altloc,
                            fullname,
                            serial_number,
                            element,
                            pqr_charge,
                            radius,
                            self.is_pqr,
                        )
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
            elif record_type == "ANISOU":
                continue
            elif record_type == "MODEL ":
                try:
                    serial_num = int(line[10:14])
                except Exception:
                    self._handle_PDB_exception(
                        "Invalid or missing model serial number", global_line_counter
                    )
                    serial_num = 0
                structure_builder.init_model(current_model_id, serial_num)
                current_model_id += 1
                model_open = 1
                current_chain_id = None
                current_residue_id = None
            elif record_type == "END   " or record_type == "CONECT":
                # End of atomic data, return the trailer
                self.line_counter += local_line_counter
                return coords_trailer[local_line_counter:]
            elif record_type == "ENDMDL":
                model_open = 0
                current_chain_id = None
                current_residue_id = None
            elif record_type == "SIGUIJ":
                continue
            elif record_type == "SIGATM":
                continue
            elif record_type not in allowed_records:
                warnings.warn(
                    f"Ignoring unrecognized record '{record_type}' at line {global_line_counter}",
                    PDBConstructionWarning,
                )
            local_line_counter += 1
        # EOF (does not end in END or CONECT)
        self.line_counter = self.line_counter + local_line_counter
        return []
