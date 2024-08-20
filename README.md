# TeV Scale Gravity Theories Branch Installation

## Steps

1. **Remove everything related to SQuIDS/nuSQuIDS**  

	Also in `/usr/local`

2. **Clone SQuIDS repository**

	```bash
	git clone https://github.com/jsalvado/SQuIDS.git
	```

3. **Clone SQuIDS 9 dimensions**

	```bash
	git clone https://github.com/pweigel/SQuIDS_dim9.git
	```

4. **Synchronize directories**

	```bash
	rsync -av ~/your/path/SQuIDS_dim9/SU_inc/ ~/your/path/SQuIDS/include/SQuIDS/SU_inc/
	```

5. **Configure SQuIDS with GSL**

	```bash
	cd SQuIDS/
	./configure --with-gsl-incdir=/usr/include --with-gsl-libdir=/usr/lib
	```

6. **Compile and test**

	```bash
	make
	make test
	sudo make install
	```

7. **Clone nuSQuIDS repo**

	```bash
	git clone -b TeVSGT https://github.com/AlanZander26/nuSQuIDS.git
	```

8. **Configure nuSQuIDS with HDF5**

	```bash
	cd nuSQuIDS/
	./configure --with-hdf5-incdir=/usr/include/hdf5/serial --with-hdf5-libdir=/usr/lib/x86_64-linux-gnu/hdf5/serial
	```

9. **Compile and test**

	```bash
	make
	make test
	make examples
	make TeVSGT
	sudo make install
	```

10. **Python binding**
Go to `Python` dir and open the Jupyter Notebook. Here you can change the model and its parameters. 

