# TeV Scale Gravity Theories Branch Installation

## Steps

1. **Remove everything related to SQuIDS/nuSQuIDS**  
	Also in `/usr/local`

2. **Clone SQuIDS repository**
	```bash
	git clone https://github.com/jsalvado/SQuIDS.git

3. **Clone SQuIDS 9 dimensions**
	```bash
	git clone https://github.com/pweigel/SQuIDS_dim9.git

4. **Synchronize directories**
	```bash
	rsync -av ~/your/path/SQuIDS_dim9/SU_inc/ ~/your/path/SQuIDS/include/SQuIDS/SU_inc/

5. **Configure SQuIDS with GSL**
	```bash
	./configure --with-gsl-incdir=/usr/include --with-gsl-libdir=/usr/lib

6. **Compile and test**
	```bash
	make
	make test
	sudo make install

7. **Clone nuSQuIDS repo**
	```bash
	git clone -b TeVSGT https://github.com/AlanZander26/nuSQuIDS.git

8. **Python binding**
Go to `Python` dir and open the Jupyter Notebook. Here you can change the model and its parameters. 

