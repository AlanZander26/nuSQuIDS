#ifndef ADD_H
#define ADD_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/tools.h>
#include <fstream>  
#include <ios>      
#include <iostream> 

namespace nusquids {

    class nuSQUIDS_ADD: public nuSQUIDS {

        private:
            unsigned int N_KK; // Number of KK states
            unsigned int dim_ADD; // Dimension of the ADD space = 3*(N_KK + 1)
            double a; // Radius of extra dimension [mu m]
            double m0; // Mass of the lightest Dirac mass [eV] (defined below m1, m2, m3. I.e, m0 = min(m1, m2, m3))
            double m1, m2, m3; // Dirac masses = eigenvalues of Dirac mass matrix
            bool NormalOrdering; // Normal neutrino mass ordering hierarchy.

            gsl_matrix_complex *W; // Transformation matrix between flavor and mass basis.

        public:

            nuSQUIDS_ADD() : nuSQUIDS(), N_KK(2), dim_ADD(3*(N_KK + 1)) {
         
                W = gsl_matrix_complex_alloc(dim_ADD, dim_ADD);
               
            }


            nuSQUIDS_ADD(marray<double,1> E_vector, unsigned int N_KK, double a, double m0, 
                bool NormalOrdering, NeutrinoType NT = both, bool iinteraction = false, 
                std::shared_ptr<CrossSectionLibrary> ncs = nullptr)
                : nuSQUIDS(E_vector, 3*(N_KK + 1), NT, iinteraction, ncs),
                    N_KK(N_KK), dim_ADD(3*(N_KK + 1)), a(a), m0(m0), NormalOrdering(NormalOrdering)
            {
                // Define standard mass splittings (in eV^2)
                double Deltaq_m21 = 7.65e-05;
                double Deltaq_m31 = 0.00247;
            
                // Set standard mixing angles
                double th01 = 0.563942, th02 = 0.154085, th12 = 0.785398;
            
                Set_MixingAngle(0, 1, th01);
                Set_MixingAngle(0, 2, th02);
                Set_MixingAngle(1, 2, th12);

                double conversion_factor = 10. / 1.98;

                std::vector<double> lambda0_list, lambda1_list, lambda2_list;

                // Solve for λ⁽⁰⁾ — lightest state
                for (unsigned int n = 0; n < N_KK+1; ++n) {
                    auto f_lambda = [=](double lambda) {
                        return lambda - M_PI * std::pow(m0 * a * conversion_factor, 2) / std::tan(M_PI * lambda);
                    };
                
                    double left = n + 1e-15;
                    double right = n + 0.5 - 1e-15;
                    double lambda_n;
                
                    try {
                        lambda_n = bisection(f_lambda, left, right);
                        lambda0_list.push_back(lambda_n);
                    } catch (const std::exception& e) {
                        std::cerr << "Error solving λ₀ for n = " << n << ": " << e.what() << std::endl;
                    }
                }
                
                // Compute zero-mode masses
                double m1_0, m2_0, m3_0;
                
                if (NormalOrdering) {
                    m1_0 = lambda0_list[0] / (a * conversion_factor);
                    m2_0 = std::sqrt(Deltaq_m21 + std::pow(m1_0, 2));
                    m3_0 = std::sqrt(Deltaq_m31 + std::pow(m1_0, 2));
                } else {
                    m3_0 = lambda0_list[0] / (a * conversion_factor);
                    m1_0 = std::sqrt(Deltaq_m31 + std::pow(m3_0, 2));
                    m2_0 = std::sqrt(Deltaq_m21 + Deltaq_m31 + std::pow(m3_0, 2));
                }
                
                // Invert to get Dirac masses
                auto invert_Dirac_mass = [&](double m0i) {
                    return std::sqrt((m0i * std::tan(M_PI * a * m0i * conversion_factor)) / (M_PI * a * conversion_factor));
                };
                                
                if (NormalOrdering) {
                    m1 = m0;
                    m2 = invert_Dirac_mass(m2_0);
                    m3 = invert_Dirac_mass(m3_0);
                } else {
                    m3 = m0;
                    m1 = invert_Dirac_mass(m1_0);
                    m2 = invert_Dirac_mass(m2_0);
                }
                
                // Now repeat transcendental root solving for m2 and m3 if NO, or m1 and m2 if IO
                double A1, A2;
                if (NormalOrdering) {
                    A1 = M_PI * std::pow(m2 * a * conversion_factor, 2);
                    A2 = M_PI * std::pow(m3 * a * conversion_factor, 2);
                } else {
                    A1 = M_PI * std::pow(m1 * a * conversion_factor, 2);
                    A2 = M_PI * std::pow(m2 * a * conversion_factor, 2);
                }
                
                // Solve λ for A1 → lambda1_list
                for (unsigned int n = 0; n < N_KK+1; ++n) {
                    auto f_lambda = [=](double lambda) {
                        return lambda - A1 / std::tan(M_PI * lambda);
                    };
                
                    double left = n + 1e-15;
                    double right = n + 0.5 - 1e-15;
                    double lambda_n;
                
                    try {
                        lambda_n = bisection(f_lambda, left, right);
                        lambda1_list.push_back(lambda_n);
                    } catch (const std::exception& e) {
                        std::cerr << "Error solving λ₁ for n = " << n << ": " << e.what() << std::endl;
                    }
                }
                
                // Solve λ for A2 → lambda2_list
                for (unsigned int n = 0; n < N_KK+1; ++n) {
                    auto f_lambda = [=](double lambda) {
                        return lambda - A2 / std::tan(M_PI * lambda);
                    };
                
                    double left = n + 1e-15;
                    double right = n + 0.5 - 1e-15;
                    double lambda_n;
                
                    try {
                        lambda_n = bisection(f_lambda, left, right);
                        lambda2_list.push_back(lambda_n);
                    } catch (const std::exception& e) {
                        std::cerr << "Error solving λ₂ for n = " << n << ": " << e.what() << std::endl;
                    }
                }             


                // Allocate mass and mixing structures
                W = gsl_matrix_complex_alloc(dim_ADD, dim_ADD);
                iniMatrices(W, th01, th02, th12);
            

                if (NormalOrdering){
                    double deltasq1 =  std::pow(lambda1_list[0]/(a * conversion_factor), 2) - std::pow(m1_0, 2);
                    Set_SquareMassDifference(1, deltasq1);
                    double deltasq2 =  std::pow(lambda2_list[0]/(a * conversion_factor), 2) - std::pow(m1_0, 2);
                    Set_SquareMassDifference(2, deltasq2);

                    // Register mass-squared differences
                    for (unsigned int n = 1; n < N_KK+1; ++n) {
                        double deltasq0 =  std::pow(lambda0_list[n]/(a * conversion_factor), 2)- std::pow(m1_0, 2);
                        Set_SquareMassDifference(3*n, deltasq0);
                        double deltasq1 =  std::pow(lambda1_list[n]/(a * conversion_factor), 2)- std::pow(m1_0, 2);
                        Set_SquareMassDifference(3*n + 1, deltasq1);
                        double deltasq2 =  std::pow(lambda2_list[n]/(a * conversion_factor), 2)- std::pow(m1_0, 2);
                        Set_SquareMassDifference(3*n + 2, deltasq2);
                    }
                } else {
                    double deltasq1 =  std::pow(lambda2_list[0]/(a * conversion_factor), 2) - std::pow(m1_0, 2);
                    Set_SquareMassDifference(1, deltasq1);
                    double deltasq2 =  std::pow(lambda0_list[0]/(a * conversion_factor), 2) - std::pow(m1_0, 2);
                    Set_SquareMassDifference(2, deltasq2);

                    // Register mass-squared differences
                    for (unsigned int n = 1; n < N_KK+1; ++n) {
                        double deltasq1 =  std::pow(lambda1_list[n]/(a * conversion_factor), 2)- std::pow(m1_0, 2);
                        Set_SquareMassDifference(3*n, deltasq1);
                        double deltasq2 =  std::pow(lambda2_list[n]/(a * conversion_factor), 2)- std::pow(m1_0, 2);
                        Set_SquareMassDifference(3*n + 1, deltasq2);
                        double deltasq0 =  std::pow(lambda0_list[n]/(a * conversion_factor), 2)- std::pow(m1_0, 2);
                        Set_SquareMassDifference(3*n + 2, deltasq0);
                    }
                }
            
               }
            

            void AddToWriteHDF5(hid_t hdf5_loc_id) const override {
            
                std::vector<double> scalars = {
                    static_cast<double>(N_KK),
                    static_cast<double>(dim_ADD),
                    a, m0, m1, m2, m3,
                    static_cast<double>(NormalOrdering)
                };
            
                hsize_t dims[1] = { scalars.size() };
                H5LTmake_dataset(hdf5_loc_id, "ADD_scalars", 1, dims, H5T_NATIVE_DOUBLE, scalars.data());
            
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

            void iniMatrices(gsl_matrix_complex*& W, double th12, double th13, double th23);

            void iniProjectors() override;

            void SetIniFlavorProyectors() override;                   
        
        
    };
} // close nusquids namespace

#endif // ADD_H
