#!/bin/bash

# Get the current directory where the script is located
TEVSGT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define the base directories based on the script's location
NU_SQUIDS_DIR="$(dirname "$TEVSGT_DIR")"

# Navigate to the 'nuSQuIDS' directory, where the makefile is
cd "$NU_SQUIDS_DIR"

# Run 'make examples_AZ'
make TeVSGT

# Navigate to the 'TeVSGT' directory
cd "$TEVSGT_DIR"

# Execute the OUT file
./tevsgt

