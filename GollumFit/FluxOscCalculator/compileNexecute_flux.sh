#!/bin/bash

# Compiler and flags
CXX=g++
CXXFLAGS="-std=c++11"

# Include directories
GSL_CFLAGS=""
GSL_LDFLAGS="-lgsl -lgslcblas -lm"
HDF5_CFLAGS="-I/home/alan/anaconda3/include"
HDF5_LDFLAGS="-L/home/alan/anaconda3/lib -lhdf5_hl -lhdf5 -L/home/alan/anaconda3/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed"
SQUIDS_CFLAGS="-Wno-abi -I/usr/local/include"
SQUIDS_LDFLAGS="-L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm"

# nuSQuIDS paths
PATH_nuSQUIDS="/home/alan/PhD/IceCube_Analysis/nuSQuIDS"
INCnuSQUIDS="$PATH_nuSQUIDS/include"
LIBnuSQUIDS="$PATH_nuSQUIDS/lib"

# Compilation and linker flags
CFLAGS="-O2 -fPIC -Ibuild -I$INCnuSQUIDS $SQUIDS_CFLAGS $GSL_CFLAGS $HDF5_CFLAGS" # Change O3 to O2
LDFLAGS="-Wl,-rpath -Wl,$LIBnuSQUIDS -L$LIBnuSQUIDS"
LDFLAGS+=" $SQUIDS_LDFLAGS $GSL_LDFLAGS $HDF5_LDFLAGS -lpthread"

# Go to nuSQuIDS directory
cd $PATH_nuSQUIDS

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
OUTPUT_PATH=$4
a=$5
m0=$6
NORMALORDERING=$7

# Compile main program with ADD.o linked
echo "Compiling the main program and linking ADD.o..."
$CXX $CXXFLAGS $CFLAGS GollumFit/FluxOscCalculator/$FILE_CPP build/ADD.o -I$PATH_nuSQUIDS/TeVSGT -lnuSQuIDS $LDFLAGS -o GollumFit/FluxOscCalculator/flux_output

if [ $? -eq 0 ]; then
    echo "Successfully compiled and linked the program."
else
    echo "Compilation failed."
    exit 1
fi

# Check if correct number of arguments are passed
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <FILE_CPP> <INPUT_FLUX> <INPUT_EARTH> <OUTPUT_PATH> <a> <m0> <NORMALORDERING>"
    exit 1
fi

# Path to the executable
EXECUTABLE="$PATH_nuSQUIDS/GollumFit/FluxOscCalculator/flux_output"

# Ensure the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found."
    exit 1
fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/alan/anaconda3/lib

# Run the executable with the provided arguments
$EXECUTABLE $INPUT_FLUX $INPUT_EARTH $OUTPUT_PATH $INDEX $a $m0 $NORMALORDERING

# Check if the program executed successfully
if [ $? -ne 0 ]; then
    echo "Error: Execution failed."
    exit 1
else
    echo "Execution completed successfully."
fi

# Erase executable
rm $EXECUTABLE
