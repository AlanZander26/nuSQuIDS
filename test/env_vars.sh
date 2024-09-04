
export CXX="g++"
export CFLAGS=" -I/home/alan/PhD/IceCube_Analysis/nuSQuIDS/inlude -Wno-abi -I/usr/local/include  -I/home/alan/anaconda3/include"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/alan/PhD/IceCube_Analysis/nuSQuIDS/lib -Wl,-rpath -Wl,/home/alan/PhD/IceCube_Analysis/nuSQuIDS/lib -lnuSQuIDS -lpthread -L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm -lgsl -lgslcblas -lm -L/home/alan/anaconda3/lib -lhdf5_hl -lhdf5 -L/home/alan/anaconda3/lib -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/alan/anaconda3/lib -Wl,-rpath-link,/home/alan/anaconda3/lib -L/home/alan/anaconda3/lib -lrt -lpthread -lz -ldl -lm -Wl,-rpath -Wl,/home/alan/anaconda3/lib"

export LD_LIBRARY_PATH="/home/alan/PhD/IceCube_Analysis/nuSQuIDS/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/home/alan/anaconda3/lib"
