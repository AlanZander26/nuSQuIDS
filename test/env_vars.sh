
export CXX="g++"
export CFLAGS=" -I/home/aurruty/Master/MasterThesis/nuSQuIDS/inlude -Wno-abi -I/usr/local/include  -I/usr/include/hdf5/serial"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/aurruty/Master/MasterThesis/nuSQuIDS/lib -Wl,-rpath -Wl,/home/aurruty/Master/MasterThesis/nuSQuIDS/lib -lnuSQuIDS -lpthread -L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm -lgsl -lgslcblas -lm -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl -lhdf5 -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial"

export LD_LIBRARY_PATH="/home/aurruty/Master/MasterThesis/nuSQuIDS/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial"
