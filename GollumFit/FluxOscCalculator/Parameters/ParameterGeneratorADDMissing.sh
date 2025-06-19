#!/bin/bash

# Path to where the .hdf5 files are stored
data_root="/data/user/azander/GolemFit/Fluxes/ADD/NO"

# Define flux components
cr_list=("GSF_1" "GSF_2" "GSF_3" "GSF_4" "GSF_5" "GSF_6")
hadron_list=("he_K+" "he_K-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_pi+" "vhe3_pi-" "vhe3_p" "vhe3_n")
fluxes=("conventional" "prompt" "astro" "${cr_list[@]}" "${hadron_list[@]}")

# Radius (a) values
a_values=(0.001000 0.001468 0.002154 0.003162 0.004642 0.006813 0.010000 0.014678 0.021544 0.031623 \
          0.046416 0.068129 0.100000 0.146780 0.215443 0.316228 0.464159 0.681292 1.000000)

# Mass (m0) values
m0_values=(0.001000 0.002154 0.004642 0.010000 0.021544 0.046416 0.100000 0.215443 0.464159 1.000000)

# Read existing .hdf5 base filenames into an array
mapfile -t existing_files < <(find "$data_root" -type f -name '*.hdf5' | sed -E 's#.*/##' | sed 's/\.hdf5$//')

# Echo how many were found
echo "Found ${#existing_files[@]} .hdf5 files in $data_root"

# Use associative array for fast lookup
declare -A existing_map
for f in "${existing_files[@]}"; do
  existing_map["$f"]=1
done

# Output file
output_file="ParametersADDMissing.txt"
> "$output_file"  # Clear file if it exists

# Generate all expected combinations and check against existing
for flux in "${fluxes[@]}"; do
  for a in "${a_values[@]}"; do
    for m0 in "${m0_values[@]}"; do
      fname="${flux}_ADD_$(printf "%.6f" "$a")_$(printf "%.6f" "$m0")"
      if [[ -z "${existing_map[$fname]}" ]]; then
        printf "%s %.6f %.6f\n" "$flux" "$a" "$m0" >> "$output_file"
      fi
    done
  done
done

echo "Missing parameters written to $output_file"
