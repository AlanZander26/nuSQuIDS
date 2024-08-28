#!/bin/bash

# Get the current directory where the script is located
TEVSGT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define the base directories based on the script's location
NU_SQUIDS_DIR="$(dirname "$TEVSGT_DIR")"

# Navigate to the 'nuSQuIDS' directory, where the makefile is
cd "$NU_SQUIDS_DIR"

# Run 'make TeVSGT'
make TeVSGT

# Navigate to the 'TeVSGT' directory
cd "$TEVSGT_DIR"

# Create a temporary copy of the 'tevsgt' binary for this job
TEMP_BINARY="./tevsgt_$$"
cp ./tevsgt "$TEMP_BINARY"

# Ensure the temporary binary is executable
chmod +x "$TEMP_BINARY"

# Run the temporary binary
"$TEMP_BINARY"

# Clean up by removing the temporary binary
rm "$TEMP_BINARY"
