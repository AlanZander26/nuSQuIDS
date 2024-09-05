#include <iostream>
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "DarkDim.h"

  
 
namespace nusquids {   

    // Returns the usual PMNS matrix.
    std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> nuSQUIDS_DarkDim::GetPMNS(double th12, double th13, double th23) {
        size_t dim_SM = 3;
        // mixing angles
        std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> th(gsl_matrix_alloc(dim_SM, dim_SM), gsl_matrix_free);  
        // cp-phases
        std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> dcp(gsl_matrix_alloc(dim_SM, dim_SM), gsl_matrix_free);

        for(unsigned int i=0; i<dim_SM; i++){
            for(unsigned int j=0; j<dim_SM; j++)
                gsl_matrix_set(th.get(),i,j,0.0);
        }
        
        for(unsigned int i=0; i<dim_SM; i++){
            for(unsigned int j=0; j<dim_SM; j++)
                gsl_matrix_set(dcp.get(),i,j,0.0);
        }

        // set mixing angles
        gsl_matrix_set(th.get(),0,1,th12);
        gsl_matrix_set(th.get(),0,2,th13);
        gsl_matrix_set(th.get(),1,2,th23);
        

        gsl_matrix_complex* U = gsl_matrix_complex_alloc(dim_SM,dim_SM);
        gsl_matrix_complex* R = gsl_matrix_complex_alloc(dim_SM,dim_SM);
        gsl_matrix_complex* dummy = gsl_matrix_complex_alloc(dim_SM,dim_SM);
        gsl_matrix_complex_set_identity(U);
        gsl_matrix_complex_set_identity(R);
        gsl_matrix_complex_set_zero(dummy);
        
        const auto unit=gsl_complex_rect(1,0);
        const auto zero=gsl_complex_rect(0,0);
        
        //construct each subspace rotation and accumulate the product
        for(size_t j=1; j<dim_SM; j++){
            for(size_t i=0; i<j; i++){
                //set up the subspace rotation
                double theta=gsl_matrix_get(th.get(),i,j);
                double delta=gsl_matrix_get(dcp.get(),i,j);
                double c=cos(theta);
                auto cp=sin(theta)*std::exp(std::complex<double>(0,-delta));
                auto cpc=-std::conj(cp);
                gsl_matrix_complex_set(R,i,i,to_gsl(c));
                gsl_matrix_complex_set(R,i,j,to_gsl(cp));
                gsl_matrix_complex_set(R,j,i,to_gsl(cpc));
                gsl_matrix_complex_set(R,j,j,to_gsl(c));
                
                //multiply this rotation onto the product from the left
                gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,unit,R,U,zero,dummy);
                std::swap(U,dummy);
                
                //clean up the rotation matrix for next iteration
                gsl_matrix_complex_set(R,i,i,unit);
                gsl_matrix_complex_set(R,i,j,zero);
                gsl_matrix_complex_set(R,j,i,zero);
                gsl_matrix_complex_set(R,j,j,unit);
            }
        }
        
        //clean up temporary matrices
        gsl_matrix_complex_free(R);
        gsl_matrix_complex_free(dummy);

        return std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>(U,gsl_matrix_complex_free);
    }



    void nuSQUIDS_DarkDim::iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th01, double th02, double th12) {
    

        std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> PMNS_temp = GetPMNS(th01, th02, th12);
        gsl_matrix_complex* PMNS = gsl_matrix_complex_alloc(3, 3);

        for(unsigned int i=0; i<3; i++){
            for(unsigned int j=0; j<3; j++){
                gsl_matrix_complex_set(PMNS, i, j, gsl_matrix_complex_get(PMNS_temp.get(),i, j));
            }
        }
       
        
        std::array<double,3> mii;
        for (int i = 0; i < 3; ++i) {
            double sqrt_factor_mii = std::sqrt((std::exp(2 * M_PI * ca[i]) - 1) / (2 * M_PI * ca[i]));
            double mi = (i == 0) ? m1 : (i == 1) ? m2 : m3;
            mii[i] = mi*sqrt_factor_mii;
        }
        

        // Compute MDc00
        gsl_matrix_complex* MDc00 = gsl_matrix_complex_alloc(3, 3);
        gsl_matrix_complex_set_zero(MDc00);

        for (int i = 0; i < 3; ++i) {
            double mi = (i == 0) ? m1 : (i == 1) ? m2 : m3;
            gsl_matrix_complex_set(MDc00, i, i, gsl_complex_rect(mi, 0));
        }

        // MDc00 = PMNS * MDc00 * PMNS^dagger
        gsl_matrix_complex* temp1 = gsl_matrix_complex_alloc(3, 3);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), PMNS, MDc00, gsl_complex_rect(0, 0), temp1);
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1, 0), temp1, PMNS, gsl_complex_rect(0, 0), MDc00);

        // aM1 construction
        gsl_matrix_complex* aM1 = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_matrix_complex_set_zero(aM1);

        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                gsl_matrix_complex_set(aM1, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(MDc00, i, j), a * (10 / 1.98)));
            }
        }

        for (int n = 1; n <= N_KK; ++n) {
            gsl_matrix_complex* MDcoff = gsl_matrix_complex_alloc(3, 3);
            gsl_matrix_complex_set_zero(MDcoff);

            for (int i = 0; i < 3; ++i) {
                double mi = mii[i];
                double sqrt_factor = std::sqrt(n * n / (n * n + ca[i] * ca[i]));
                gsl_matrix_complex_set(MDcoff, i, i, gsl_complex_rect(mi * sqrt_factor, 0));
            }

            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), PMNS, MDcoff, gsl_complex_rect(0, 0), temp1);
            gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1, 0), temp1, PMNS, gsl_complex_rect(0, 0), MDcoff);

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    gsl_matrix_complex_set(aM1, 3 * n + i, j, gsl_complex_mul_real(gsl_matrix_complex_get(MDcoff, i, j), std::sqrt(2) * a * (10 / 1.98)));
                }
            }

            gsl_matrix_complex_free(MDcoff);
        }
    

        // aM2 construction
        gsl_matrix_complex* aM2 = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_matrix_complex_set_zero(aM2);

        gsl_matrix_complex* aMD2 = gsl_matrix_complex_alloc(3, 3);

        for (int n = 1; n <= N_KK; ++n) {
            
            gsl_matrix_complex_set_zero(aMD2);

            for (int i = 0; i < 3; ++i) {
                double sqrt_factor = std::sqrt(n * n + ca[i] * ca[i]);
                gsl_matrix_complex_set(aMD2, i, i, gsl_complex_rect(sqrt_factor, 0));
            }

            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), PMNS, aMD2, gsl_complex_rect(0, 0), temp1);
            gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1, 0), temp1, PMNS, gsl_complex_rect(0, 0), aMD2);

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    gsl_matrix_complex_set(aM2, 3 * n + i, 3 * n + j, gsl_matrix_complex_get(aMD2, i, j));
                }
            }
        }


        gsl_matrix_complex_free(aMD2);
        gsl_matrix_complex_free(temp1);
 

        // Combine aM1 and aM2 to form aM
        gsl_matrix_complex* aM = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_matrix_complex_memcpy(aM, aM1);
        gsl_matrix_complex_add(aM, aM2);


        // Calculate aaMM = aM^dagger * aM
        gsl_matrix_complex* aaMM = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1, 0), aM, aM, gsl_complex_rect(0, 0), aaMM);

        // Cleanup
        gsl_matrix_complex_free(PMNS);
        gsl_matrix_complex_free(MDc00);
        gsl_matrix_complex_free(aM1);
        gsl_matrix_complex_free(aM2);
        gsl_matrix_complex_free(aM);


        gsl_matrix_complex *M2 = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1)); // Mass matrix squared in eV^2.
        gsl_complex value, newValue;
        for (int i = 0; i < 3 * (N_KK + 1); ++i) {
            for (int j = 0; j < 3 * (N_KK + 1); ++j) {
                value = gsl_matrix_complex_get(aaMM, i, j);
                newValue = gsl_complex_div_real(value, std::pow(a, 2)/3.92e-2); 
                gsl_matrix_complex_set(M2, i, j, newValue);
            }
        }    

        gsl_matrix_complex_free(aaMM);


        gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(3*(N_KK + 1));
        gsl_vector* Lambdaq_temp = gsl_vector_alloc(3*(N_KK + 1));
        
        // Perform eigenvalue decomposition
        gsl_eigen_hermv(M2, Lambdaq_temp, W, w);

        


        // Free workspace
        gsl_eigen_hermv_free(w);

        // Copy eigenvalues to a vector and sort it
    std::vector<double> Lambdaq_vec(Lambdaq_temp->size);
    std::vector<size_t> index_map(Lambdaq_temp->size);
    for (size_t i = 0; i < Lambdaq_temp->size; ++i) {
        Lambdaq_vec[i] = gsl_vector_get(Lambdaq_temp, i);
        index_map[i] = i; // Initialize index map
    }
  

    // Sort Lambdaq_vec and maintain index mapping
    std::sort(index_map.begin(), index_map.end(), [&Lambdaq_vec](size_t i1, size_t i2) {
        return Lambdaq_vec[i1] < Lambdaq_vec[i2];
    });


    // Create sorted matrix W_sorted
    gsl_matrix_complex *W_sorted = gsl_matrix_complex_alloc(W->size1, W->size2);

    // Rearrange columns of W according to sorted eigenvalues and save to W_sorted
    for (size_t j = 0; j < index_map.size(); ++j) {
        size_t sorted_index = index_map[j];
        for (size_t i = 0; i < W->size1; ++i) {
            gsl_complex value = gsl_matrix_complex_get(W, i, sorted_index);
            gsl_matrix_complex_set(W_sorted, i, j, value);
        }
    }

    // Copy the sorted values back to Lambdaq_temp
    for (size_t i = 0; i < Lambdaq_temp->size; ++i) {
        gsl_vector_set(Lambdaq_temp, i, Lambdaq_vec[index_map[i]]);
    }

    // Calculate differences and store in Lambdaq
    for (size_t i = 0; i < Lambdaq_temp->size - 1; ++i) {
        double value = gsl_vector_get(Lambdaq_temp, i + 1) - gsl_vector_get(Lambdaq_temp, 0);
        gsl_vector_set(Lambdaq, i, value);
    }
 
        // Create a temporary matrix to store the sorted columns of W

        // Copy the sorted W back to the original W
        gsl_matrix_complex_memcpy(W, W_sorted);

        gsl_vector_free(Lambdaq_temp);
        gsl_matrix_complex_free(W_sorted);
   

    }

/*
 void nuSQUIDS_DarkDim::iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th01, double th02, double th12) {
    

        std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> PMNS_temp = GetPMNS(th01, th02, th12);
        gsl_matrix_complex* PMNS = gsl_matrix_complex_alloc(3, 3);

        for(unsigned int i=0; i<3; i++){
            for(unsigned int j=0; j<3; j++){
                gsl_matrix_complex_set(PMNS, i, j, gsl_matrix_complex_get(PMNS_temp.get(),i, j));
            }
        }
       
        
        std::array<double,3> mii;
        for (int i = 0; i < 3; ++i) {
            double sqrt_factor_mii = std::sqrt((std::exp(2 * M_PI * ca[i]) - 1) / (2 * M_PI * ca[i]));
            double mi = (i == 0) ? m1 : (i == 1) ? m2 : m3;
            mii[i] = mi*sqrt_factor_mii;
        }
        

        // Compute MDc00
        gsl_matrix_complex* MDc00_diag = gsl_matrix_complex_alloc(3, 3);
        gsl_matrix_complex_set_zero(MDc00_diag);

        for (int i = 0; i < 3; ++i) {
            double mi = (i == 0) ? m1 : (i == 1) ? m2 : m3;
            gsl_matrix_complex_set(MDc00_diag, i, i, gsl_complex_rect(mi, 0));
        }

        // MDc00 = PMNS * MDc00_diag * PMNS^dagger
        gsl_matrix_complex* MDc00 = gsl_matrix_complex_alloc(3, 3);
        gsl_matrix_complex_set_zero(MDc00);
        gsl_matrix_complex* temp1 = gsl_matrix_complex_alloc(3, 3);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), PMNS, MDc00_diag, gsl_complex_rect(0, 0), temp1);
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1, 0), temp1, PMNS, gsl_complex_rect(0, 0), MDc00);
        gsl_matrix_complex_free(temp1);
        gsl_matrix_complex_free(MDc00_diag);

        // aM1 construction
        gsl_matrix_complex* aM1 = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_matrix_complex_set_zero(aM1);

        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                gsl_matrix_complex_set(aM1, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(MDc00, i, j), a * (10 / 1.98)));
            }
        }

        gsl_matrix_complex* temp2 = gsl_matrix_complex_alloc(3, 3);
        for (int n = 1; n <= N_KK; ++n) {
            gsl_matrix_complex* MDcoff = gsl_matrix_complex_alloc(3, 3);
            gsl_matrix_complex_set_zero(MDcoff);
            

            for (int i = 0; i < 3; ++i) {
                double mi = mii[i];
                double sqrt_factor = std::sqrt(n * n / (n * n + ca[i] * ca[i]));
                gsl_matrix_complex_set(MDcoff, i, i, gsl_complex_rect(mi * sqrt_factor, 0));
            }

            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), PMNS, MDcoff, gsl_complex_rect(0, 0), temp2);
            gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1, 0), temp2, PMNS, gsl_complex_rect(0, 0), MDcoff);

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    gsl_matrix_complex_set(aM1, 3 * n + i, j, gsl_complex_mul_real(gsl_matrix_complex_get(MDcoff, i, j), std::sqrt(2) * a * (10 / 1.98)));
                }
            }

            gsl_matrix_complex_free(MDcoff);
        }
        gsl_matrix_complex_free(temp2);

        // aM2 construction
        gsl_matrix_complex* aM2 = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_matrix_complex_set_zero(aM2);

        gsl_matrix_complex* aMD2 = gsl_matrix_complex_alloc(3, 3);
        gsl_matrix_complex* temp3 = gsl_matrix_complex_alloc(3, 3);  

        for (int n = 1; n <= N_KK; ++n) {
            gsl_matrix_complex_set_zero(aMD2);

            for (int i = 0; i < 3; ++i) {
                double sqrt_factor = std::sqrt(n * n + ca[i] * ca[i]);
                gsl_matrix_complex_set(aMD2, i, i, gsl_complex_rect(sqrt_factor, 0));
            }

            // Multiply PMNS * aMD2
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), PMNS, aMD2, gsl_complex_rect(0, 0), temp3);

            // Multiply the result with the conjugate transpose of PMNS
            gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1, 0), temp3, PMNS, gsl_complex_rect(0, 0), aMD2);

            // Set the corresponding block in aM2
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    gsl_matrix_complex_set(aM2, 3 * n + i, 3 * n + j, gsl_matrix_complex_get(aMD2, i, j));
                }
            }
        }

// Free allocated matrices
gsl_matrix_complex_free(aMD2);
gsl_matrix_complex_free(temp3);

        // Combine aM1 and aM2 to form aM
        gsl_matrix_complex* aM = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_matrix_complex_memcpy(aM, aM1);
        gsl_matrix_complex_add(aM, aM2);

        // Calculate aaMM = aM^dagger * aM
        gsl_matrix_complex* aaMM = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1));
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1, 0), aM, aM, gsl_complex_rect(0, 0), aaMM);

        // Cleanup
        gsl_matrix_complex_free(PMNS);
        gsl_matrix_complex_free(MDc00);
        gsl_matrix_complex_free(aM1);
        gsl_matrix_complex_free(aM2);
        gsl_matrix_complex_free(aM);


        gsl_matrix_complex *M2 = gsl_matrix_complex_alloc(3 * (N_KK + 1), 3 * (N_KK + 1)); // Mass matrix squared in eV^2.
        gsl_complex value, newValue;
        for (int i = 0; i < 3 * (N_KK + 1); ++i) {
            for (int j = 0; j < 3 * (N_KK + 1); ++j) {
                value = gsl_matrix_complex_get(aaMM, i, j);
                newValue = gsl_complex_div_real(value, std::pow(a, 2)/3.92e-2); 
                gsl_matrix_complex_set(M2, i, j, newValue);
            }
        }    

        gsl_matrix_complex_free(aaMM);


        gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(3*(N_KK + 1));
        gsl_vector* Lambdaq_temp = gsl_vector_alloc(3*(N_KK + 1));
        
        // Perform eigenvalue decomposition
        gsl_eigen_hermv(M2, Lambdaq_temp, W, w);


        // Free workspace
        gsl_eigen_hermv_free(w);

        // Copy eigenvalues to a vector and sort it
    std::vector<double> Lambdaq_vec(Lambdaq_temp->size);
    std::vector<size_t> index_map(Lambdaq_temp->size);
    for (size_t i = 0; i < Lambdaq_temp->size; ++i) {
        Lambdaq_vec[i] = gsl_vector_get(Lambdaq_temp, i);
        index_map[i] = i; // Initialize index map
    }
  

    // Sort Lambdaq_vec and maintain index mapping
    std::sort(index_map.begin(), index_map.end(), [&Lambdaq_vec](size_t i1, size_t i2) {
        return Lambdaq_vec[i1] < Lambdaq_vec[i2];
    });


    // Create sorted matrix W_sorted
    gsl_matrix_complex *W_sorted = gsl_matrix_complex_alloc(W->size1, W->size2);

    // Rearrange columns of W according to sorted eigenvalues and save to W_sorted
    for (size_t j = 0; j < index_map.size(); ++j) {
        size_t sorted_index = index_map[j];
        for (size_t i = 0; i < W->size1; ++i) {
            gsl_complex value = gsl_matrix_complex_get(W, i, sorted_index);
            gsl_matrix_complex_set(W_sorted, i, j, value);
        }
    }

    // Copy the sorted values back to Lambdaq_temp
    for (size_t i = 0; i < Lambdaq_temp->size; ++i) {
        gsl_vector_set(Lambdaq_temp, i, Lambdaq_vec[index_map[i]]);
    }

    // Calculate differences and store in Lambdaq
    for (size_t i = 0; i < Lambdaq_temp->size - 1; ++i) {
        double value = gsl_vector_get(Lambdaq_temp, i + 1) - gsl_vector_get(Lambdaq_temp, 0);
        gsl_vector_set(Lambdaq, i, value);
    }
      
        // Create a temporary matrix to store the sorted columns of W

        // Copy the sorted W back to the original W
        gsl_matrix_complex_memcpy(W, W_sorted);

        gsl_vector_free(Lambdaq_temp);
        gsl_matrix_complex_free(W_sorted);
   

    }
    */



    void nuSQUIDS_DarkDim::iniProjectors(){

        b0_proj.resize(std::vector<size_t>{numneu});
        for(unsigned int flv = 0; flv < numneu; flv++){
        b0_proj[flv] = squids::SU_vector::Projector(nsun,flv);
        }

        b1_proj.resize(std::vector<size_t>{nrhos,numneu});
        for(unsigned int rho = 0; rho < nrhos; rho++){
            for(unsigned int flv = 0; flv < numneu; flv++){
                b1_proj[rho][flv] = squids::SU_vector::Projector(nsun,flv);

                AntineutrinoCPFix(rho);
                // Rotate with transformation matrix
                b1_proj[rho][flv] = b1_proj[rho][flv].Rotate(W);
                AntineutrinoCPFix(rho);
            }
        }

        evol_b0_proj.resize(std::vector<size_t>{nrhos,numneu,ne});
        evol_b1_proj.resize(std::vector<size_t>{nrhos,numneu,ne});
        for(unsigned int rho = 0; rho < nrhos; rho++){
            for(unsigned int flv = 0; flv < numneu; flv++){
                for(unsigned int e1 = 0; e1 < ne; e1++){
                evol_b0_proj[rho][flv][e1] = squids::SU_vector::Projector(nsun,flv);
                evol_b1_proj[rho][flv][e1] = squids::SU_vector::Projector(nsun,flv);

                AntineutrinoCPFix(rho);
                // Rotate with transformation matrix
                evol_b1_proj[rho][flv][e1] = evol_b1_proj[rho][flv][e1].Rotate(W);
                AntineutrinoCPFix(rho);
                }
            }
        }
    } 


  void nuSQUIDS_DarkDim::SetIniFlavorProyectors(){

    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        for(unsigned int e1 = 0; e1 < ne; e1++){
          evol_b1_proj[rho][flv][e1] = b0_proj[flv];

          AntineutrinoCPFix(rho);
          evol_b1_proj[rho][flv][e1] = evol_b1_proj[rho][flv][e1].Rotate(W);
          AntineutrinoCPFix(rho);
        }
        b1_proj[rho][flv] = b0_proj[flv];

        AntineutrinoCPFix(rho);
        b1_proj[rho][flv] = b1_proj[rho][flv].Rotate(W);
        AntineutrinoCPFix(rho);
      }
    }
  }
/*
  std::array<double, 3> nuSQUIDS_DarkDim::calculate_lambda() {
    M_Pl = 2.435e18 // Reduced Planck mass [GeV]
    std::array<double, 3> lambda_i;
    std::array<double, 3> m = {m1, m2, m3};

    for (size_t i = 0; i < 3; ++i) {
        double exponent = 2 * M_PI * ca[i];
        double numerator = exp(exponent) - 1;
        double denominator = exponent;
        lambda_i[i] = m[i] * M_Pl / Mf * sqrt(numerator / denominator);
    }

    return lambda_i;
}*/

}

