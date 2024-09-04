
export CXX="g++"
<<<<<<< HEAD
export CFLAGS=" -I/home/alan/PhD/IceCube_Analysis/nuSQuIDS/inlude -Wno-abi -I/usr/local/include  -I/home/alan/anaconda3/include"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/alan/PhD/IceCube_Analysis/nuSQuIDS/lib -Wl,-rpath -Wl,/home/alan/PhD/IceCube_Analysis/nuSQuIDS/lib -lnuSQuIDS -lpthread -L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm -lgsl -lgslcblas -lm -L/home/alan/anaconda3/lib -lhdf5_hl -lhdf5 -L/home/alan/anaconda3/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/alan/anaconda3/lib -Wl,-rpath-link,/home/alan/anaconda3/lib -L/home/alan/anaconda3/lib -lrt -lpthread -lz -ldl -lm -Wl,-rpath -Wl,/home/alan/anaconda3/lib"

export LD_LIBRARY_PATH="/home/alan/PhD/IceCube_Analysis/nuSQuIDS/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/home/alan/anaconda3/lib"
=======
export CFLAGS=" -I/home/agurruty/Projects/MasterThesis/nuSQuIDS/inlude -Wno-abi -I/usr/local/include  -I/usr/include/hdf5/serial"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/agurruty/Projects/MasterThesis/nuSQuIDS/lib -Wl,-rpath -Wl,/home/agurruty/Projects/MasterThesis/nuSQuIDS/lib -lnuSQuIDS -lpthread -L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm -lgsl -lgslcblas -lm -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl"

export LD_LIBRARY_PATH="/home/agurruty/Projects/MasterThesis/nuSQuIDS/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu/hdf5/serial:/home/agurruty/Projects/SciComputing/root_install/lib"
>>>>>>> d75229f7f4a6e21fe7524df41cdf3d61b0449dd6
