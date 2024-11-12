#!/bin/bash

MODEL=$1
NORMALORDERING=$2
shift 2 # Shift to make $@ contain only PARAMS

PARAMS="$@" # Capture all remaining arguments as PARAMS

# Define an array with the different flux types
flux_types=("conventional" "prompt" "astro")

# Loop over each flux type and execute the script with the required arguments
for flux in "${flux_types[@]}"; do
    echo "Running compileNexecute_flux.sh with FLUX_TYPE: $flux"
    bash compileNexecute_flux.sh "$flux" "$NORMALORDERING" "$MODEL" $PARAMS
    echo "Finished running for FLUX_TYPE: $flux"
    echo
done

#Example: bash generate_fluxes_point.sh "ADD" true 0.500000 0.000000
