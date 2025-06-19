#!/bin/bash

USER_COBALT=$1
NORMALORDERING=$2
MODEL=$3

shift 3 # Shift to make $@ contain only PARAMS

PARAMS="$@" # Capture all remaining arguments as PARAMS

# Define the GOLEMDIR environment variable
GOLEMDIR=$(realpath ../../../../)

# Source the environment setup
source ${GOLEMDIR}/local/setup.sh

FLUXDIR=/data/user/${USER_COBALT}/GolemFit

# Define the hadron and CR lists
cr_list=("GSF_1" "GSF_2" "GSF_3" "GSF_4" "GSF_5" "GSF_6")
hadron_list=("he_K+" "he_K-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_pi+" "vhe3_pi-" "vhe3_p" "vhe3_n")
# Complete list:
#("he_K+" "he_K-" "he_n" "he_p" "he_pi+" "he_pi-" "le_K+" "le_K-" \
#        "le_pi+" "le_pi-" "vhe1_pi+" "vhe1_pi-" "vhe3_K+" "vhe3_K-" "vhe3_n" \
#       "vhe3_p" "vhe3_pi+" "vhe3_pi-")

# Combine all flux types into a single array
flux_types=("conventional" "prompt" "astro" "${cr_list[@]}" "${hadron_list[@]}")

# Loop over each flux type and execute the script with the required arguments
for flux in "${flux_types[@]}"; do
    echo "Running compileNexecute_flux.sh with FLUX_TYPE: $flux"
    bash compileNexecute_flux.sh "$FLUXDIR" "$flux" "$NORMALORDERING" "$MODEL" $PARAMS
    echo "Finished running for FLUX_TYPE: $flux"
    echo
done

#Example: bash generate_fluxes_point_cobalt.sh "azander" true "ADD" 0.500000 0.000000
