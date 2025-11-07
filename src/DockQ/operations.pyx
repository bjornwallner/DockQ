#cython: language_level=3
import numpy as np
cimport numpy as np
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
def residue_distances(float [:,:] atom_coordinates1, float [:,:] atom_coordinates2, long [:] atoms_per_res1, long [:] atoms_per_res2):
    cdef:
        int i, j, x, y, i_atoms, j_atoms, cum_i_atoms, cum_j_atoms, cum_i_atoms_end, cum_j_atoms_end
        int n_res_i = atoms_per_res1.shape[0]
        int n_res_j = atoms_per_res2.shape[0]
        float this_d, min_d
        float [:,:] res_distances = np.zeros((n_res_i, n_res_j), dtype=np.float32)

    cum_i_atoms = 0
    for i in range(n_res_i):
        i_atoms = atoms_per_res1[i]
        cum_i_atoms_end = cum_i_atoms + i_atoms
        cum_j_atoms = 0
        for j in range(n_res_j):
            j_atoms = atoms_per_res2[j]
            min_d = 100000.0
            cum_j_atoms_end = cum_j_atoms + j_atoms
            for x in range(cum_i_atoms, cum_i_atoms_end):
                for y in range(cum_j_atoms, cum_j_atoms_end):
                    this_d = (atom_coordinates1[x][0] - atom_coordinates2[y][0])**2 + (atom_coordinates1[x][1] - atom_coordinates2[y][1])**2 + (atom_coordinates1[x][2] - atom_coordinates2[y][2])**2
                    if this_d < min_d:
                        min_d = this_d
                        if min_d > 1000.0:
                            break
                if min_d > 1000.0:
                    break
            res_distances[i, j] = min_d
            cum_j_atoms = cum_j_atoms + j_atoms
        cum_i_atoms = cum_i_atoms + i_atoms

    return res_distances

@cython.boundscheck(False)
@cython.wraparound(False)
def get_fnat_stats(float [:,:] model_res_distances, float [:,:] native_res_distances, float threshold=5.0):
    cdef:
        int native_shape_0 = native_res_distances.shape[0]
        int native_shape_1 = native_res_distances.shape[1]
        int i, j
        int n_native_contacts = 0
        int n_model_contacts = 0
        int n_shared_contacts = 0
        int n_non_native_contacts = 0
        float threshold_squared

    threshold_squared = threshold * threshold
    for i in range(native_shape_0):
        for j in range(native_shape_1):
            if native_res_distances[i, j] < threshold_squared:
                n_native_contacts += 1
                if model_res_distances[i, j] < threshold_squared:
                    n_shared_contacts += 1
            if model_res_distances[i, j] < threshold_squared:
                n_model_contacts += 1
                if native_res_distances[i, j] >= threshold_squared:
                    n_non_native_contacts += 1

    return (
        n_shared_contacts,
        n_non_native_contacts,
        n_native_contacts,
        n_model_contacts,
    )
