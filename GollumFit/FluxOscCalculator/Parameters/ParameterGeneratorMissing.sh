#!/bin/bash

# Define the base directory
BASE_DIR="/data/user/azander/GolemFit/Fluxes/SM_Copies/IO"

# Define the output file
output_file="Parameters_missing.txt"
> "$output_file"  # Clear the file if it exists

# Define the hadron and CR lists
cr_list=("GSF_1" "GSF_2" "GSF_3" "GSF_4" "GSF_5" "GSF_6")
hadron_list=("he_K+" "he_K-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_pi+" "vhe3_pi-" "vhe3_p" "vhe3_n")

# Loop through all directories in the base directory
for dir in "$BASE_DIR"/*/; do
    if [[ -d "$dir" ]]; then
        # Extract parameters from the directory name correctly
        dir_name=$(basename "$dir")
        read -r p1 p2 p3 <<< $(echo "$dir_name" | awk -F'_' '{printf "%.6f %.6f %.6f", $3, $4, $5}')

  # Check for missing subdirectories and log if they are absent
for subdir in "conventional" "prompt" "astro"; do
    if [[ ! -d "$dir/$subdir" ]]; then
        if [[ "$p2" != "499.999995" ]]; then
            echo "$subdir $p1 $p2 $p3" >> "$output_file"
        fi
    else
        # Count only .hdf5 files
        hdf5_count=$(find "$dir/$subdir" -maxdepth 1 -type f -name "*.hdf5" | wc -l)
        if [[ "$hdf5_count" -ne 1 ]]; then
            if [[ "$p2" != "499.999995" ]]; then
                echo "$subdir $p1 $p2 $p3" >> "$output_file"
            fi
        fi
    fi
done

        # Check errors directory for 16 or 24 .hdf5 files
        error_count=$(find "$dir/errors" -maxdepth 1 -type f -name "*.hdf5" | wc -l)
        if [[ "$error_count" -ne 16 && "$error_count" -ne 24 ]]; then
            for flux in "${cr_list[@]}" "${hadron_list[@]}"; do
                expected_file="$dir/errors/${flux}_SM_Copies_${p1}_${p2}_${p3}.hdf5"
                if [[ ! -f "$expected_file" ]]; then
                    if [[ "$p2" != "499.999995" ]]; then
                        echo "$flux $p1 $p2 $p3" >> "$output_file"
                    fi
                fi
            done
        fi
        
    fi
done

echo "Check complete. Missing parameters logged in $output_file."




