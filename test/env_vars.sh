
export CXX="g++"
export CFLAGS=" -I/home/agurruty/Projects/MasterThesis/nuSQuIDS/inlude -Wno-abi -I/usr/local/include  -I/usr/include/hdf5/serial"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/agurruty/Projects/MasterThesis/nuSQuIDS/lib -Wl,-rpath -Wl,/home/agurruty/Projects/MasterThesis/nuSQuIDS/lib -lnuSQuIDS -lpthread -L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm -lgsl -lgslcblas -lm -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl"

export LD_LIBRARY_PATH="/home/agurruty/Projects/MasterThesis/nuSQuIDS/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial:/home/agurruty/Projects/SciComputing/root_install/lib"
