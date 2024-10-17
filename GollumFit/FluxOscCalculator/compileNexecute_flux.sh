#!/bin/bash

export LD_LIBRARY_PATH=~/GOLEMSOURCE/local/lib:$LD_LIBRARY_PATH
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

# Compile ADD.cpp into object file
echo "Compiling ADD.cpp into ADD.o..."
$CXX $CXXFLAGS $CFLAGS -c TeVSGT/ADD/ADD.cpp -o build/ADD.o

if [ $? -eq 0 ]; then
    echo "Successfully compiled ADD.cpp."
else
    echo "Compilation of ADD.cpp failed."
    exit 1
fi

# Assign arguments to variables
FILE_CPP=$1
INPUT_FLUX=$2
INPUT_EARTH=$3
FLUX_NAME=$4
a=$5
m0=$6
NORMALORDERING=$7

point="ADD_${a}_${m0}"
OUTPUT_PATH=${SCRIPT_DIR}/Fluxes/$point

mkdir -p $OUTPUT_PATH

NAME_EXECUTABLE=flux_output_${FILE_CPP}_${point}

# Compile main program with ADD.o linked
echo "Compiling the main program and linking ADD.o..."
$CXX $CXXFLAGS $CFLAGS GollumFit/FluxOscCalculator/$FILE_CPP build/ADD.o -I$PATH_nuSQUIDS/TeVSGT -lnuSQuIDS $LDFLAGS -o GollumFit/FluxOscCalculator/$NAME_EXECUTABLE

if [ $? -eq 0 ]; then
    echo "Successfully compiled and linked the program."
else
    echo "Compilation failed."
    exit 1
fi

# Check if correct number of arguments are passed
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <FILE_CPP> <INPUT_FLUX> <INPUT_EARTH> <FLUX_NAME> <a> <m0> <NORMALORDERING>"
    exit 1
fi

# Path to the executable
EXECUTABLE="$PATH_nuSQUIDS/GollumFit/FluxOscCalculator/$NAME_EXECUTABLE"

# Ensure the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found."
    exit 1
fi


# Run the executable with the provided arguments
$EXECUTABLE $INPUT_FLUX $INPUT_EARTH $OUTPUT_PATH $FLUX_NAME $a $m0 $NORMALORDERING

# Check if the program executed successfully
if [ $? -ne 0 ]; then
    echo "Error: Execution failed."
    exit 1
else
    echo "Execution completed successfully."
fi

# Erase executable
rm $EXECUTABLE

# Example of usage: bash compileNexecute_flux.sh prompt_atmospheric_flux.cpp GollumFit/FluxOscCalculator/AIRS_mceq121_pr_flux_2011_sib_HG_E3.dat GollumFit/FluxOscCalculator/EARTH_MODEL_PREM.dat "new_ddm" 0.500000 0.000000 true
