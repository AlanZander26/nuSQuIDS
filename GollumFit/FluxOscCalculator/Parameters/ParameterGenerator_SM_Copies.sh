#!/bin/bash

m0=0.000000

# Define the functions for logspace
logspace() {
    local start=$1
    local end=$2
    local points=$3
    awk -v start="$start" -v end="$end" -v points="$points" 'BEGIN {
        for (i = 0; i < points; i++) {
            value = 10 ^ (start + i * (end - start) / (points - 1))
            printf "%.6f\n", value
        }
    }'
}

# Define the hadron and CR lists
cr_list=("GSF_1" "GSF_2" "GSF_3" "GSF_4" "GSF_5" "GSF_6")
hadron_list=("he_K+" "he_K-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_pi+" "vhe3_pi-" "vhe3_p" "vhe3_n")
fluxes=("astro" "prompt" "conventional" "${cr_list[@]}" "${hadron_list[@]}")

# Generate N_values and mu_values
N_values=$(logspace 1 $(awk 'BEGIN {print log(500)/log(10)}') 20)
mu_values=$(logspace 0 $(awk 'BEGIN {print log(500)/log(10)}') 20)

# Create the param list
> nusquidsParameters.txt
for flux in "${fluxes[@]}"; do
    for N in $N_values; do
        for mu in $mu_values; do
            # Skip mu values >= 499.99
            if (( $(echo "$mu == 499.999995" | bc -l) )); then
                continue
            fi
            echo "$flux $N $mu $m0" >> nusquidsParameters.txt
        done
    done
done


