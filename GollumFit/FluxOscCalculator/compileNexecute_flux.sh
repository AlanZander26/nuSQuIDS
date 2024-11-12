#!/bin/bash

# Assign common arguments to variables
FLUX_TYPE=$1 # Flux type ("conventional", "prompt" or "astro")
INPUT_FLUX=$2
INPUT_EARTH=$3
NORMALORDERING=$4
MODEL=$5  # Model type
shift 5   # Shift past the first 5 arguments

Model=${MODEL^^}
sed -i "1s/.*/#define USE_${Model}/" calculate_flux.cpp

# Parse remaining arguments based on MODEL
case "$MODEL" in
  "ADD")
    # ADD model: expects `a` and `m0`
    if [ "$#" -ne 2 ]; then
      echo "Error: ADD model requires exactly 2 additional arguments: a and m0."
      exit 1
    fi
    a=$1
    m0=$2
    ;;

  "SM_Copies")
    # SM_Copies model: expects `N`, `mu`, and `m0`
    if [ "$#" -ne 3 ]; then
      echo "Error: SM_Copies model requires exactly 3 additional arguments: N, mu, and m0."
      exit 1
    fi
    N=$1
    mu=$2
    m0=$3
    ;;

    "SM")
    # SM: expects no more argument
    if [ "$#" -ne 0 ]; then
      echo "Error: SM requires no additional argument."
      exit 1
    fi
    ;;

  *)
    echo "Error: Unsupported model type '$MODEL'."
    exit 1
    ;;
esac

# Display parsed variables (for debugging purposes)
echo "FLUX_TYPE: $FLUX_TYPE"
echo "INPUT_FLUX: $INPUT_FLUX"
echo "INPUT_EARTH: $INPUT_EARTH"
echo "NORMALORDERING: $NORMALORDERING"
echo "MODEL: $MODEL"

# Display model-specific variables
if [ "$MODEL" == "ADD" ]; then
    echo "ADD Model Parameters:"
    echo "  a: $a"
    echo "  m0: $m0"
elif [ "$MODEL" == "SM_Copies" ]; then
    echo "SM_Copies Model Parameters:"
    echo "  N: $N"
    echo "  mu: $mu"
    echo "  m0: $m0"
fi


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

if [ "$MODEL" != "SM" ]; then
# Compile .cpp file into object file
echo "Compiling $MODEL.cpp into $MODEL.o..."
$CXX $CXXFLAGS $CFLAGS -c TeVSGT/$MODEL/$MODEL.cpp -o build/$MODEL.o
fi


if [ $? -eq 0 ] && [ "$MODEL" != "SM" ]; then
    echo "Successfully compiled $MODEL.cpp."
elif [ "$MODEL" == "SM" ]; then
    :
else
    echo "Compilation of $MODEL.cpp failed."
    exit 1
fi

# Define model-dependent parameters
if [ "$MODEL" == "ADD" ]; then
    point="${MODEL}_${a}_${m0}"
elif [ "$MODEL" == "SM_Copies" ]; then
    point="${MODEL}_${N}_${mu}_${m0}"
elif [ "$MODEL" == "SM" ]; then
    points="${MODEL}"  # No additional parameters for the Standard Model
else
    echo "Error: Unknown model $MODEL"
    exit 1
fi

OUTPUT_PATH=${SCRIPT_DIR}/Fluxes/$MODEL/$point

mkdir -p $OUTPUT_PATH

NAME_EXECUTABLE=flux_output_${FLUX_TYPE}_${point}

# Compile main program.
if [ "$MODEL" == "SM" ]; then
    FILEO=""
    echo "Compiling the main program..."
else
    FILEO="build/$MODEL.o"
    echo "Compiling the main program and linking $MODEL.o..."
fi
$CXX $CXXFLAGS $CFLAGS GollumFit/FluxOscCalculator/calculate_flux.cpp $FILEO -I$PATH_nuSQUIDS/TeVSGT -lnuSQuIDS $LDFLAGS -o GollumFit/FluxOscCalculator/$NAME_EXECUTABLE

if [ $? -eq 0 ]; then
    echo "Successfully compiled and linked the program."
else
    echo "Compilation failed."
    exit 1
fi

# Path to the executable
EXECUTABLE="$PATH_nuSQUIDS/GollumFit/FluxOscCalculator/$NAME_EXECUTABLE"

# Ensure the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE not found."
    exit 1
fi

# Define model-dependent parameters
if [ "$MODEL" == "ADD" ]; then
    PARAMS="$a $m0"
elif [ "$MODEL" == "SM_Copies" ]; then
    PARAMS="$N $mu $m0"
elif [ "$MODEL" == "SM" ]; then
    PARAMS=""  # No additional parameters for the Standard Model
else
    echo "Error: Unknown model $MODEL"
    exit 1
fi

# Run the executable with the provided arguments
$EXECUTABLE $FLUX_TYPE $INPUT_FLUX $INPUT_EARTH $OUTPUT_PATH $NORMALORDERING $PARAMS

# Check if the program executed successfully
if [ $? -ne 0 ]; then
    echo "Error: Execution failed."
    exit 1
else
    echo "Execution completed successfully."
fi

# Erase executable
rm $EXECUTABLE

# Example of usage: bash compileNexecute_flux.sh "conventional" GollumFit/FluxOscCalculator/Data/v0.6.0_nodeis/ddm_conv_bestfit.dat GollumFit/FluxOscCalculator/Data/EARTH_MODEL_PREM.dat true "ADD" 0.500000 0.000000 
