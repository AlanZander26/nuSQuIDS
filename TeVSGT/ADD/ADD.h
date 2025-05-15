#ifndef ADD_H
#define ADD_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids {

    class nuSQUIDS_ADD: public nuSQUIDS {

        private:
            unsigned int N_KK; // Number of KK states
            unsigned int dim_ADD; // Dimension of the ADD space = 3*(N_KK + 1)
            double a; // Radius of extra dimension [mu m]
            double m0; // Mass of the lightest neutrino state [eV]
            double m1, m2, m3; // Masses of light neutrinos.
            bool NormalOrdering; // Normal neutrino mass ordering hierarchy.
            double th01, th02, th12;

    

            gsl_vector *Lambdaq; // Vector with "squared masses" (lambdas, see paper Machado et.al.) in ascending order.
            gsl_matrix_complex *W; // Transformation matrix between flavor and mass basis.

        public:

            nuSQUIDS_ADD(marray<double,1> E_vector, unsigned int N_KK, double a, double m0, 
            bool NormalOrdering, NeutrinoType NT = both, bool iinteraction = false, 
            std::shared_ptr<CrossSectionLibrary> ncs = nullptr) : N_KK(N_KK), dim_ADD(3*(N_KK + 1)), a(a), m0(m0), NormalOrdering(NormalOrdering), 
            nuSQUIDS(E_vector, 3*(N_KK + 1), NT, iinteraction, ncs)
            {
                //===============================
                // set mixing angles and squared mass differences   //
                //===============================        
                
                double Deltaq_m21 = 7.65e-05; // Change the squared mass difference here in eV^2
                double Deltaq_m31 = 0.00247; // Change the squared mass difference here in eV^2          
                

                if (NormalOrdering) {
                    m1 = m0;
                    m2 = std::sqrt(Deltaq_m21 + std::pow(m0, 2));
                    m3 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
                } else {
                    m1 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
                    m2 = std::sqrt(Deltaq_m21 + Deltaq_m31 + std::pow(m0, 2));
                    m3 = m0;
                }   

                //Lambdaq = gsl_vector_alloc(dim_ADD-1);
                //W = gsl_matrix_complex_alloc(dim_ADD, dim_ADD);
                //iniMatrices(Lambdaq, W, th01, th02, th12);
                iniProjectors();

                Set_MixingAngle(0,1,th01);
                Set_MixingAngle(0,2,th02);
                Set_MixingAngle(1,2,th12);

                for (int j = 0; j < dim_ADD-1; j++) {
                    Set_SquareMassDifference(j+1,gsl_vector_get(Lambdaq, j));
                }
            }

 


/*
        void AddToWriteHDF5(hid_t hdf5_loc_id) const {
            // Step 1: Get dimensions of the W matrix
            unsigned int rows = W->size1; // Number of rows in W
            unsigned int cols = W->size2; // Number of columns in W
            hsize_t W_dim[2] = {rows, cols * 2}; // HDF5 stores real and imaginary parts separately

            // Step 2: Flatten the W matrix into a vector of doubles
            std::vector<double> W_flat(rows * cols * 2); // Each complex number has 2 doubles (real, imag)

            for (unsigned int i = 0; i < rows; ++i) {
                for (unsigned int j = 0; j < cols; ++j) {
                    gsl_complex z = gsl_matrix_complex_get(W, i, j);
                    W_flat[(i * cols + j) * 2]     = GSL_REAL(z); // Real part
                    W_flat[(i * cols + j) * 2 + 1] = GSL_IMAG(z); // Imaginary part
                }
            }

            // Step 3: Write the flattened matrix to the HDF5 file
            H5LTmake_dataset(hdf5_loc_id, "W_matrix", 2, W_dim, H5T_NATIVE_DOUBLE, W_flat.data());
        }


        void AddToReadHDF5(hid_t hdf5_loc_id) {
            // Step 1: Get dataset dimensions
            hsize_t dims[2]; // Array to hold dimensions [rows, cols*2]
            H5LTget_dataset_info(hdf5_loc_id, "W_matrix", dims, nullptr, nullptr);

            unsigned int rows = dims[0];       // Number of rows
            unsigned int cols = dims[1] / 2;   // Number of columns (divide by 2 because of real and imag parts)

            // Step 2: Read flattened W matrix from the HDF5 file
            std::unique_ptr<double[]> W_data(new double[dims[0] * dims[1]]);
            H5LTread_dataset_double(hdf5_loc_id, "W_matrix", W_data.get());

            // Step 3: Allocate and reconstruct the gsl_matrix_complex W
            W = gsl_matrix_complex_alloc(rows, cols); // Ensure W is properly allocated

            for (unsigned int i = 0; i < rows; ++i) {
                for (unsigned int j = 0; j < cols; ++j) {
                    double real_part = W_data[(i * cols + j) * 2];       // Real part
                    double imag_part = W_data[(i * cols + j) * 2 + 1];   // Imaginary part
                    gsl_complex z = gsl_complex_rect(real_part, imag_part); // Construct complex number
                    gsl_matrix_complex_set(W, i, j, z);                 // Set value in W
                }
            }
        }
*/
            std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetPMNS(double th12 = 0.563942, double th13 = 0.154085, double th23 = 0.785398);

            void iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th12, double th13, double th23);

            void iniProjectors() override;

            void SetIniFlavorProyectors() override;      
        
        
    };
} // close nusquids namespace

#endif // ADD_H