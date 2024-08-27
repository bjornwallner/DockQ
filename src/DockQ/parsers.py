import warnings
import numpy as np
import Bio
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBParser import as_handle
from Bio.SeqUtils import seq1

custom_map = {"MSE": "M", "CME": "C"}


class MMCIFParser(Bio.PDB.MMCIFParser):
    def get_structure(
        self,
        structure_id,
        filename,
        chains=[],
        parse_hetatms=False,
        auth_chains=True,
        model_number=0,
    ):
        """Return the structure.

        Arguments:
         - structure_id - string, the id that will be used for the structure
         - filename - name of mmCIF file, OR an open text mode file handle

        """
        self.auth_chains = auth_chains
        self.auth_residues = True
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            self._mmcif_dict = MMCIF2Dict(filename)
            sequences, is_het = self._build_structure(
                structure_id, chains, parse_hetatms=parse_hetatms
            )
            self._structure_builder.set_header(self._get_header())

        structure = self._structure_builder.get_structure()
        model = structure[model_number]
        for chain in model:
            chain.sequence = sequences[chain.id]
            chain.is_het = is_het[chain.id]
        return model

    def _build_structure(self, structure_id, chains, parse_hetatms):
        # two special chars as placeholders in the mmCIF format
        # for item values that cannot be explicitly assigned
        # see: pdbx/mmcif syntax web page
        _unassigned = {".", "?"}

        mmcif_dict = self._mmcif_dict
        sequences = {}
        is_het = {}
        atom_serial_list = mmcif_dict["_atom_site.id"]
        atom_id_list = mmcif_dict["_atom_site.label_atom_id"]
        residue_id_list = mmcif_dict["_atom_site.label_comp_id"]
        try:
            element_list = mmcif_dict["_atom_site.type_symbol"]
        except KeyError:
            element_list = None
        if self.auth_chains:
            chain_id_list = mmcif_dict["_atom_site.auth_asym_id"]
        else:
            chain_id_list = mmcif_dict["_atom_site.label_asym_id"]
        x_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_x"]]
        y_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_y"]]
        z_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_z"]]
        alt_list = mmcif_dict["_atom_site.label_alt_id"]
        icode_list = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
        b_factor_list = mmcif_dict["_atom_site.B_iso_or_equiv"]
        occupancy_list = mmcif_dict["_atom_site.occupancy"]
        fieldname_list = mmcif_dict["_atom_site.group_PDB"]
        try:
            serial_list = [int(n) for n in mmcif_dict["_atom_site.pdbx_PDB_model_num"]]
        except KeyError:
            # No model number column
            serial_list = None
        except ValueError:
            # Invalid model number (malformed file)
            raise PDBConstructionException("Invalid model number") from None
        try:
            aniso_u11 = mmcif_dict["_atom_site_anisotrop.U[1][1]"]
            aniso_u12 = mmcif_dict["_atom_site_anisotrop.U[1][2]"]
            aniso_u13 = mmcif_dict["_atom_site_anisotrop.U[1][3]"]
            aniso_u22 = mmcif_dict["_atom_site_anisotrop.U[2][2]"]
            aniso_u23 = mmcif_dict["_atom_site_anisotrop.U[2][3]"]
            aniso_u33 = mmcif_dict["_atom_site_anisotrop.U[3][3]"]
            aniso_flag = 1
        except KeyError:
            # no anisotropic B factors
            aniso_flag = 0

        if self.auth_residues:
            # if auth_seq_id is present, we use this.
            # Otherwise label_seq_id is used.
            if "_atom_site.auth_seq_id" in mmcif_dict:
                seq_id_list = mmcif_dict["_atom_site.auth_seq_id"]
            else:
                seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        else:
            seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        # Now loop over atoms and build the structure
        current_chain_id = None
        current_residue_id = None
        current_resname = None
        structure_builder = self._structure_builder
        structure_builder.init_structure(structure_id)
        structure_builder.init_seg(" ")
        # Historically, Biopython PDB parser uses model_id to mean array index
        # so serial_id means the Model ID specified in the file
        current_model_id = -1
        current_serial_id = -1
        for i in range(len(atom_id_list)):
            chainid = chain_id_list[i]
            if chains and chainid not in chains:
                continue
            fieldname = fieldname_list[i]
            if fieldname == "HETATM" and not parse_hetatms:
                continue
            element = element_list[i].upper() if element_list else None
            if element == "H":
                continue
            # set the line_counter for 'ATOM' lines only and not
            # as a global line counter found in the PDBParser()
            structure_builder.set_line_counter(i)

            # Try coercing serial to int, for compatibility with PDBParser
            # But do not quit if it fails. mmCIF format specs allow strings.
            try:
                serial = int(atom_serial_list[i])
            except ValueError:
                serial = atom_serial_list[i]
                warnings.warn(
                    "PDBConstructionWarning: Some atom serial numbers are not numerical",
                    PDBConstructionWarning,
                )

            x = x_list[i]
            y = y_list[i]
            z = z_list[i]
            resname = residue_id_list[i]
            altloc = alt_list[i]
            if altloc in _unassigned:
                altloc = " "
            resseq = seq_id_list[i]
            if resseq == ".":
                # Non-existing residue ID
                try:
                    msg_resseq = mmcif_dict["_atom_site.auth_seq_id"][i]
                    msg = f"Non-existing residue ID in chain '{chainid}', residue '{msg_resseq}'"
                except (KeyError, IndexError):
                    msg = f"Non-existing residue ID in chain '{chainid}'"
                warnings.warn(
                    "PDBConstructionWarning: " + msg,
                    PDBConstructionWarning,
                )
                continue
            int_resseq = int(resseq)
            icode = icode_list[i]
            if icode in _unassigned:
                icode = " "
            name = atom_id_list[i]
            # occupancy & B factor
            try:
                tempfactor = float(b_factor_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing B factor") from None
            try:
                occupancy = float(occupancy_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing occupancy") from None
            hetatm_flag = " " if fieldname != "HETATM" else "H"

            resseq = (hetatm_flag, int_resseq, icode)

            if serial_list is not None:
                # model column exists; use it
                serial_id = serial_list[i]
                if current_serial_id != serial_id:
                    # if serial changes, update it and start new model
                    current_serial_id = serial_id
                    current_model_id += 1
                    structure_builder.init_model(current_model_id, current_serial_id)
                    current_chain_id = None
                    current_residue_id = None
                    current_resname = None
            else:
                # no explicit model column; initialize single model
                structure_builder.init_model(current_model_id)

            if current_chain_id != chainid:
                current_chain_id = chainid
                if current_chain_id not in sequences:
                    sequences[current_chain_id] = ""
                if current_chain_id not in is_het:
                    is_het[current_chain_id] = None
                structure_builder.init_chain(current_chain_id)
                current_residue_id = None
                current_resname = None

            if current_residue_id != resseq or current_resname != resname:
                current_residue_id = resseq
                current_resname = resname
                if hetatm_flag == " ":
                    resname1 = (
                        seq1(current_resname, custom_map=custom_map)
                        if len(current_resname) == 3
                        else current_resname[-1]
                        if (len(current_resname) == 2)
                        else current_resname
                    )
                    sequences[current_chain_id] += resname1
                else:
                    sequences[current_chain_id] = resname
                    is_het[current_chain_id] = resname
                structure_builder.init_residue(resname, hetatm_flag, int_resseq, icode)

            coord = np.array((x, y, z), "f")

            structure_builder.init_atom(
                name,
                coord,
                tempfactor,
                occupancy,
                altloc,
                name,
                serial_number=serial,
                element=element,
            )
            if aniso_flag == 1 and i < len(aniso_u11):
                u = (
                    aniso_u11[i],
                    aniso_u12[i],
                    aniso_u13[i],
                    aniso_u22[i],
                    aniso_u23[i],
                    aniso_u33[i],
                )
                mapped_anisou = [float(_) for _ in u]
                anisou_array = np.array(mapped_anisou, "f")
                structure_builder.set_anisou(anisou_array)
        # Now try to set the cell

        try:
            a = float(mmcif_dict["_cell.length_a"][0])
            b = float(mmcif_dict["_cell.length_b"][0])
            c = float(mmcif_dict["_cell.length_c"][0])
            alpha = float(mmcif_dict["_cell.angle_alpha"][0])
            beta = float(mmcif_dict["_cell.angle_beta"][0])
            gamma = float(mmcif_dict["_cell.angle_gamma"][0])
            cell = np.array((a, b, c, alpha, beta, gamma), "f")
            spacegroup = mmcif_dict["_symmetry.space_group_name_H-M"][0]
            spacegroup = spacegroup[1:-1]  # get rid of quotes!!
            if spacegroup is None:
                raise Exception
            structure_builder.set_symmetry(spacegroup, cell)
        except Exception:
            pass  # no cell found, so just ignore
        return sequences, is_het


class PDBParser(Bio.PDB.PDBParser):
    def get_structure(self, id, file, chains, parse_hetatms, model_number=0):
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
                sequences, is_het = self._parse(
                    lines, chains, parse_hetatms=parse_hetatms
                )
            self.structure_builder.set_header(self.header)
            # Return the Structure instance
            structure = self.structure_builder.get_structure()
            model = structure[model_number]

            for chain in model:
                chain.sequence = sequences[chain.id]
                chain.is_het = is_het[chain.id]
        return model

    def _parse(self, header_coords_trailer, chains, parse_hetatms):
        """Parse the PDB file (PRIVATE)."""
        # Extract the header; return the rest of the file
        self.header, coords_trailer = self._get_header(header_coords_trailer)
        # Parse the atomic data; return the PDB file trailer
        self.trailer, sequences, is_het = self._parse_coordinates(
            coords_trailer, chains, parse_hetatms
        )
        return sequences, is_het

    def _parse_coordinates(self, coords_trailer, chains=[], parse_hetatms=False):
        """Parse the atomic data in the PDB file (PRIVATE)."""
        allowed_records = {
            "ATOM  ",
            "HETATM",
            "MODEL ",
            "ENDMDL",
            "TER   ",
        }
        sequences = {}
        is_het = {}
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
            elif record_type == "HETATM" and not parse_hetatms:
                continue
            elif record_type == "ATOM  " or record_type == "HETATM":
                # Initialize the Model - there was no explicit MODEL record
                if not model_open:
                    structure_builder.init_model(current_model_id)
                    current_model_id += 1
                    model_open = 1
                chainid = line[21]
                if chains and chainid not in chains:
                    continue
                element = line[76:78].strip().upper()
                if element == "H":
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
                hetatm_flag = " "
                if record_type == "HETATM":  # hetero atom flag
                    # if a small molecule and the name matches what we're looking for
                    hetatm_flag = "H"

                try:
                    serial_number = int(line[6:11])
                except Exception:
                    serial_number = 0
                resseq = int(line[22:26].split()[0])  # sequence identifier
                icode = line[26]  # insertion code
                residue_id = (hetatm_flag, resseq, icode)
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

                segid = line[72:76]
                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id = residue_id
                    current_resname = resname
                    if current_chain_id not in sequences:
                        sequences[current_chain_id] = ""
                    if current_chain_id not in is_het:
                        is_het[current_chain_id] = None
                    try:
                        structure_builder.init_residue(
                            resname, hetatm_flag, resseq, icode
                        )
                        if hetatm_flag == " ":
                            resname1 = (
                                seq1(current_resname, custom_map=custom_map)
                                if len(current_resname) == 3
                                else current_resname[-1]
                                if (len(current_resname) == 2)
                                else current_resname
                            )
                            sequences[current_chain_id] = resname1
                        else:
                            sequences[current_chain_id] = current_resname
                            is_het[current_chain_id] = current_resname
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetatm_flag, resseq, icode
                        )
                        if hetatm_flag == " ":
                            resname1 = (
                                seq1(current_resname, custom_map=custom_map)
                                if len(current_resname) == 3
                                else current_resname[-1]
                                if (len(current_resname) == 2)
                                else current_resname
                            )
                            sequences[current_chain_id] += resname1
                        else:
                            sequences[current_chain_id] = current_resname
                            is_het[current_chain_id] = current_resname
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)

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
                return coords_trailer[local_line_counter:], sequences, is_het
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
        return [], sequences, is_het
