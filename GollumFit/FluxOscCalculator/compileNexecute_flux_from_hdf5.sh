#!/bin/bash

# Check that the user passed exactly one argument
if [ "$#" -ne 1 ]; then
    echo "Usage: bash compileNexecute_flux_from_hdf5.sh /path/to/file.hdf5"
    exit 1
fi

# HDF5 file path (full)
HDF5_FILE="$1"

export LD_LIBRARY_PATH=~/ADD_GOLEMSOURCE/local/lib:$LD_LIBRARY_PATH
# Compiler and flags
CXX=g++
CXXFLAGS="-std=c++11"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set PATH_nuSQUIDS to two levels up from the script's directory
PATH_nuSQUIDS=$(realpath "$SCRIPT_DIR/../..")

# Go to nuSQuIDS directory
cd $PATH_nuSQUIDS

eval $(sed -n '/^GSL_CFLAGS=/,/^LIBnuSQUIDS=/p' $PATH_nuSQUIDS/Makefile \
    | sed 's/=\(.*\)/="\1"/' \
    | sed 's/\$(\([a-zA-Z_][a-zA-Z0-9_]*\))/\${\1}/g'
)

# Compilation and linker flags
CFLAGS="-O2 -fPIC -Ibuild -I$INCnuSQUIDS $SQUIDS_CFLAGS $GSL_CFLAGS $HDF5_CFLAGS" # Change O3 to O2
LDFLAGS="-Wl,-rpath -Wl,$LIBnuSQUIDS -L$LIBnuSQUIDS"
LDFLAGS+=" $SQUIDS_LDFLAGS $GSL_LDFLAGS $HDF5_LDFLAGS -lpthread"

if [ "$MODEL" != "SM" ]; then
# Compile .cpp file into object file
echo "Compiling SM_Copies.cpp into SM_Copies.o..." # $MODEL.cpp into $MODEL.o..."
$CXX $CXXFLAGS $CFLAGS -c TeVSGT/SM_Copies/SM_Copies.cpp -o build/SM_Copies.o # TeVSGT/$MODEL/$MODEL.cpp -o build/$MODEL.o
fi

# Output binary name
EXEC="flux_reader"

FILEO="build/SM_Copies.o"

# Compile the program
$CXX $CXXFLAGS $CFLAGS GollumFit/FluxOscCalculator/compute_flux_hdf5.cpp $FILEO -I$PATH_nuSQUIDS/TeVSGT -lnuSQuIDS $LDFLAGS -o GollumFit/FluxOscCalculator/$EXEC

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed."
    exit 1
fi

echo "Compilation succeeded."

# Run the compiled program with the HDF5 file
EXECUTABLE="$PATH_nuSQUIDS/GollumFit/FluxOscCalculator/$EXEC"

# Ensure the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found."
    exit 1
fi

# Run the executable with the provided arguments
$EXECUTABLE $HDF5_FILE > "$PATH_nuSQUIDS/GollumFit/FluxOscCalculator/flux_from_hdf5.log" 

# Check if the program executed successfully
if [ $? -ne 0 ]; then
    echo "Error: Execution failed."
    exit 1
else
    echo "Execution completed successfully."
fi

# Erase executable
rm $EXECUTABLE

echo "Done."
