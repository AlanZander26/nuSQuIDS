#!/bin/bash

export PATH=$PATH:~/GOLEMSOURCE/local/bin
source ~/GOLEMSOURCE/local/setup.sh

# Read the first non-comment line from fixedParameters.txt
fixed_line=$(grep -v '^#' Parameters/fixedParameters.txt | head -n 1)

# Split the line into an array of parameters
read -r USER_COBALT NORMALORDERING MODEL <<< "$fixed_line"

FLUXDIR=/data/user/${USER_COBALT}/GolemFit

flux=$1

shift 1 # Shift to make $@ contain only PARAMS

PARAMS="$@" # Capture all remaining arguments as PARAMS

bash compileNexecute_flux.sh "$FLUXDIR" "$flux" "$NORMALORDERING" "$MODEL" $PARAMS
