#cython: language_level=3
import numpy as np
cimport numpy as np
import cython
from libc.math cimport sqrt


def residue_distances(float [:,:] atom_coordinates, long [:] atoms_per_res1, long [:] atoms_per_res2, float threshold=5.0):
    cdef:
        int i, j, x, y
        int n_res_i = len(atoms_per_res1)
        int n_res_j = len(atoms_per_res2)
        float min_d = 100.0
        float this_d
        float [:,:] res_distances = np.zeros((n_res_i, n_res_j), dtype=np.float32)


    for i in range(n_res_i):
        n_atoms_res_i = atoms_per_res1[i]
        for j in range(n_res_j):
            n_atoms_res_j = atoms_per_res2[j]
            
            for x in range(n_atoms_res_i):
                for y in range(n_atoms_res_j):
                    this_d = sqrt((atom_coordinates[x][0] - atom_coordinates[y][0])**2 + (atom_coordinates[x][1] - atom_coordinates[y][1])**2 + (atom_coordinates[x][2] - atom_coordinates[y][2])**2)
                    if this_d < min_d:
                        min_d = this_d
                if min_d < threshold:
                    break
            res_distances[i, j] = min_d
    
    return res_distances


def get_fnat_stats(float [:,:] model_res_distances, float [:,:] native_res_distances, float threshold=5.0):
    cdef:
        int model_shape_0 = model_res_distances.shape[0]
        int model_shape_1 = model_res_distances.shape[1]
        int native_shape_0 = native_res_distances.shape[0]
        int native_shape_1 = native_res_distances.shape[1]
        int i, j
        int n_native_contacts = 0
        int n_model_contacts = 0
        int n_shared_contacts = 0
        int n_non_native_contacts = 0
        #unsigned char [:,:] model_contacts = np.zeros((model_shape_0, model_shape_1), dtype=np.uint8)
        #unsigned char [:,:] native_contacts = np.zeros((native_shape_0, native_shape_1), dtype=np.uint8)
    
    for i in range(native_shape_0):
        for j in range(native_shape_1):
            if native_res_distances[i, j] < threshold:
                n_native_contacts += 1
                if model_res_distances[i, j] < threshold:
                    n_shared_contacts += 1
            if model_res_distances[i, j] < threshold:
                n_model_contacts += 1
                if native_res_distances[i, j] >= threshold:
                    n_non_native_contacts += 1

    return (
        n_shared_contacts,
        n_non_native_contacts,
        n_native_contacts,
        n_model_contacts,
    )
