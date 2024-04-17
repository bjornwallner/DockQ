import numpy as np


def get_distances_across_chains(model_A_atoms, model_B_atoms):

    distances = ((model_A_atoms[:, None] - model_B_atoms[None, :]) ** 2).sum(-1)

    return distances


def atom_distances_to_residue_distances(atom_distances, atoms_per_res1, atoms_per_res2):
    res_distances = np.zeros((len(atoms_per_res1), len(atoms_per_res2)))

    cum_i_atoms = 0
    for i, i_atoms in enumerate(atoms_per_res1):
        cum_j_atoms = 0
        for j, j_atoms in enumerate(atoms_per_res2):
            res_distances[i, j] = atom_distances[
                cum_i_atoms : cum_i_atoms + i_atoms, cum_j_atoms : cum_j_atoms + j_atoms
            ].min()
            cum_j_atoms += j_atoms
        cum_i_atoms += i_atoms
    return res_distances


def residue_distances(
    atom_coordinates1, atom_coordinates2, atoms_per_res1, atoms_per_res2
):
    atom_distances = get_distances_across_chains(atom_coordinates1, atom_coordinates2)
    res_distances = atom_distances_to_residue_distances(
        atom_distances, atoms_per_res1, atoms_per_res2
    )

    return res_distances


def get_fnat_stats(model_res_distances, native_res_distances, threshold=5.0):
    native_contacts = native_res_distances < threshold ** 2
    model_contacts = model_res_distances < threshold ** 2
    n_native_contacts = np.sum(native_contacts)
    n_model_contacts = np.sum(model_contacts)
    n_shared_contacts = np.sum(model_contacts * native_contacts)
    n_non_native_contacts = np.sum(model_contacts * (1 - native_contacts))
    return (
        n_shared_contacts,
        n_non_native_contacts,
        n_native_contacts,
        n_model_contacts,
    )
