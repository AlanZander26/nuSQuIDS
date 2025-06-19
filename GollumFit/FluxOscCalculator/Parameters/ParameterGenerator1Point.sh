#!/bin/bash

# Check if at least one argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 val1 [val2 val3 ...]"
    exit 1
fi

# Define flux lists
cr_list=("GSF_1" "GSF_2" "GSF_3" "GSF_4" "GSF_5" "GSF_6")
hadron_list=("he_K+" "he_K-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_pi+" "vhe3_pi-" "vhe3_p" "vhe3_n")
flux_types=("conventional" "prompt" "astro" "${cr_list[@]}" "${hadron_list[@]}")

# Output file
output_file="Parameters1Point.txt"
> "$output_file"  # clear or create the output file

# Loop over all flux types
for flux in "${flux_types[@]}"; do
    echo "$flux $*" >> "$output_file"
done

echo "Output written to $output_file"
