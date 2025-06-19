#!/bin/bash

# Define CRs and Hadrons
cr_list=("GSF_1" "GSF_2" "GSF_3" "GSF_4" "GSF_5" "GSF_6")
hadron_list=("he_K+" "he_K-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_pi+" "vhe3_pi-" "vhe3_p" "vhe3_n")
fluxes=("conventional" "prompt" "astro" "${cr_list[@]}" "${hadron_list[@]}") for the 1st run "astro" always fails

# Define 'a' values (radius) with 6 decimal places
a_values=(0.001000 0.001468 0.002154 0.003162 0.004642 0.006813 0.010000 0.014678 0.021544 0.031623 \
          0.046416 0.068129 0.100000 0.146780 0.215443 0.316228 0.464159 0.681292 1.000000)

# Define 'm0' values with 6 decimal places
m0_values=(0.001000 0.002154 0.004642 0.010000 0.021544 0.046416 0.100000 0.215443 0.464159 1.000000)

# Output file
output_file="Parameters_ADD.txt"
> "$output_file"  # Clear existing file if it exists

# Generate combinations
for flux in "${fluxes[@]}"; do
  for a in "${a_values[@]}"; do
    for m0 in "${m0_values[@]}"; do
      printf "%s %.6f %.6f\n" "$flux" "$a" "$m0" >> "$output_file"
    done
  done
done

echo "Output written to $output_file"
