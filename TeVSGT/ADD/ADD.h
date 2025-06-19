#ifndef ADD_H
#define ADD_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <fstream>  
#include <ios>      
#include <iostream> 

namespace nusquids {

    class nuSQUIDS_ADD: public nuSQUIDS {

        private:
            unsigned int N_KK; // Number of KK states
            unsigned int dim_ADD; // Dimension of the ADD space = 3*(N_KK + 1)
            double a; // Radius of extra dimension [mu m]
            double m0; // Mass of the lightest neutrino state [eV]
            double m1, m2, m3; // Masses of light neutrinos.
            bool NormalOrdering; // Normal neutrino mass ordering hierarchy.

            gsl_vector *Lambdaq; // Vector with "squared masses" (lambdas, see paper Machado et.al.) in ascending order.
            gsl_matrix_complex *W; // Transformation matrix between flavor and mass basis.

        public:

            nuSQUIDS_ADD() : nuSQUIDS(), N_KK(2), dim_ADD(3*(N_KK + 1)) {
         
                Lambdaq = gsl_vector_alloc(dim_ADD-1);
                W = gsl_matrix_complex_alloc(dim_ADD, dim_ADD);
               
            }

            nuSQUIDS_ADD(marray<double,1> E_vector, unsigned int N_KK, double a, double m0, 
            bool NormalOrdering, NeutrinoType NT = both, bool iinteraction = false, 
            std::shared_ptr<CrossSectionLibrary> ncs = nullptr) : nuSQUIDS(E_vector, 3*(N_KK + 1), NT, iinteraction, ncs), N_KK(N_KK), dim_ADD(3*(N_KK + 1)), a(a), m0(m0), NormalOrdering(NormalOrdering)
            {
                //===============================
                // set mixing angles and squared mass differences   //
                //===============================        
                
                double Deltaq_m21 = 7.65e-05; // Change the squared mass difference here in eV^2
                double Deltaq_m31 = 0.00247; // Change the squared mass difference here in eV^2
                double th01=0.563942, th02=0.154085, th12=0.785398;
            
                Set_MixingAngle(0,1,th01);
                Set_MixingAngle(0,2,th02);
                Set_MixingAngle(1,2,th12);

                if (NormalOrdering) {
                    m1 = m0;
                    m2 = std::sqrt(Deltaq_m21 + std::pow(m0, 2));
                    m3 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
                } else {
                    m1 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
                    m2 = std::sqrt(Deltaq_m21 + Deltaq_m31 + std::pow(m0, 2));
                    m3 = m0;
                }   

                Lambdaq = gsl_vector_alloc(dim_ADD-1);
                W = gsl_matrix_complex_alloc(dim_ADD, dim_ADD);
                iniMatrices(Lambdaq, W, th01, th02, th12);


                for (unsigned int j = 0; j < dim_ADD-1; j++) {
                    Set_SquareMassDifference(j+1,gsl_vector_get(Lambdaq, j));
                }
            }


            void AddToWriteHDF5(hid_t hdf5_loc_id) const override {
                std::cout << "[nuSQUIDS_ADD] Writing scalar attributes as dataset to HDF5..." << std::endl;
            
                std::vector<double> scalars = {
                    static_cast<double>(N_KK),
                    static_cast<double>(dim_ADD),
                    a, m0, m1, m2, m3,
                    static_cast<double>(NormalOrdering)
                };
            
                hsize_t dims[1] = { scalars.size() };
                H5LTmake_dataset(hdf5_loc_id, "ADD_scalars", 1, dims, H5T_NATIVE_DOUBLE, scalars.data());
            
                std::cout << "  Wrote dataset 'ADD_scalars' with values:" << std::endl;
                for (size_t i = 0; i < scalars.size(); ++i) {
                    std::cout << "    [" << i << "] = " << scalars[i] << std::endl;
                }
            
                // Lambdaq
                if (Lambdaq != nullptr) {
                    hsize_t dims[1] = { Lambdaq->size };
                    H5LTmake_dataset(hdf5_loc_id, "Lambdaq", 1, dims, H5T_NATIVE_DOUBLE, Lambdaq->data);
                }
            
                // W matrix
                if (W != nullptr) {
                    hsize_t dims[2] = { W->size1, W->size2 };
                    std::vector<double> W_real(W->size1 * W->size2);
                    std::vector<double> W_imag(W->size1 * W->size2);
            
                    for (size_t i = 0; i < W->size1; ++i) {
                        for (size_t j = 0; j < W->size2; ++j) {
                            gsl_complex z = gsl_matrix_complex_get(W, i, j);
                            W_real[i * W->size2 + j] = GSL_REAL(z);
                            W_imag[i * W->size2 + j] = GSL_IMAG(z);
                        }
                    }
            
                    H5LTmake_dataset(hdf5_loc_id, "W_real", 2, dims, H5T_NATIVE_DOUBLE, W_real.data());
                    H5LTmake_dataset(hdf5_loc_id, "W_imag", 2, dims, H5T_NATIVE_DOUBLE, W_imag.data());
                }
            }
            

            void AddToReadHDF5(hid_t hdf5_loc_id) override {
                std::cout << "[nuSQUIDS_ADD] AddToReadHDF5() called." << std::endl;

                // Read scalar values from dataset
                hsize_t dims[1];
                H5LTget_dataset_info(hdf5_loc_id, "ADD_scalars", dims, nullptr, nullptr);
            
                if (dims[0] < 8) {
                    throw std::runtime_error("ADD_scalars dataset is too small");
                }
            
                std::vector<double> scalars(dims[0]);
                H5LTread_dataset_double(hdf5_loc_id, "ADD_scalars", scalars.data());
            
                N_KK = static_cast<unsigned int>(scalars[0]);
                dim_ADD = static_cast<unsigned int>(scalars[1]);
                a = scalars[2];
                m0 = scalars[3];
                m1 = scalars[4];
                m2 = scalars[5];
                m3 = scalars[6];
                NormalOrdering = (scalars[7] > 0.5);
            
                std::cout << "[nuSQUIDS_ADD] Read dataset 'ADD_scalars' with values:" << std::endl;
                for (size_t i = 0; i < scalars.size(); ++i) {
                    std::cout << "    [" << i << "] = " << scalars[i] << std::endl;
                }
            
                // Lambdaq
                H5LTget_dataset_info(hdf5_loc_id, "Lambdaq", dims, nullptr, nullptr);
                if (Lambdaq) gsl_vector_free(Lambdaq);
                Lambdaq = gsl_vector_alloc(dims[0]);
                H5LTread_dataset_double(hdf5_loc_id, "Lambdaq", Lambdaq->data);
            
                // W
                hsize_t mdims[2];
                H5LTget_dataset_info(hdf5_loc_id, "W_real", mdims, nullptr, nullptr);
                if (W) gsl_matrix_complex_free(W);
                W = gsl_matrix_complex_alloc(mdims[0], mdims[1]);
            
                std::vector<double> W_real(mdims[0] * mdims[1]);
                std::vector<double> W_imag(mdims[0] * mdims[1]);
            
                H5LTread_dataset_double(hdf5_loc_id, "W_real", W_real.data());
                H5LTread_dataset_double(hdf5_loc_id, "W_imag", W_imag.data());
            
                for (size_t i = 0; i < mdims[0]; ++i) {
                    for (size_t j = 0; j < mdims[1]; ++j) {
                        gsl_complex z;
                        GSL_SET_COMPLEX(&z, W_real[i * mdims[1] + j], W_imag[i * mdims[1] + j]);
                        gsl_matrix_complex_set(W, i, j, z);
                    }
                }
            }
            

            std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetPMNS(double th12 = 0.563942, double th13 = 0.154085, double th23 = 0.785398);

            void iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th12, double th13, double th23);

            void iniProjectors() override;

            void SetIniFlavorProyectors() override;                   
        
        
    };
} // close nusquids namespace

#endif // ADD_H
