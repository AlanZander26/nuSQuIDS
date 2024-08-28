#!/bin/bash

# Define the ranges over which to iterate
model_types=("ADD")
n_mesh_values=(10)
n_final_mesh_values=(10)
E_min_values=(1)
E_max_values=(1e5)

# Create the params_ranges_list dynamically
params_ranges_list=()

# Generate combinations for the parameter ranges
for i in $(seq 0.1 0.1 0.3); do
    for j in $(seq 0 0.05 0.05); do
        params_ranges_list+=("[[$i,$j]]")
    done
done

# GNU Parallel for distributing tasks to multiple CPUs/nodes
parallel -j 6 --ungroup --link \
    'python3 generateFluxes.py {1} {2} {3} {4} {5} {6}' ::: \
    "${model_types[@]}" ::: \
    "${params_ranges_list[@]}" ::: \
    "${n_mesh_values[@]}" ::: \
    "${n_final_mesh_values[@]}" ::: \
    "${E_min_values[@]}" ::: \
    "${E_max_values[@]}"
