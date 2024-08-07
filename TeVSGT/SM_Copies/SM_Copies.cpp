#include <complex>
#include <gsl/gsl_blas.h> 

#include "SM_Copies.h"

namespace nusquids {
  

  auto to_gsl=[](const std::complex<double>& c)->gsl_complex {
    return(gsl_complex_rect(c.real(),c.imag()));
  };

  // Returns the usual PMNS matrix.
  std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> nuSQUIDS_SM_Copies::GetPMNS(double th12, double th13, double th23) {
    // mixing angles
    std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> th(gsl_matrix_alloc(dim, dim), gsl_matrix_free);  
    // cp-phases
    std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> dcp(gsl_matrix_alloc(dim, dim), gsl_matrix_free);

    for(unsigned int i=0; i<dim; i++){
          for(unsigned int j=0; j<dim; j++)
              gsl_matrix_set(th.get(),i,j,0.0);
      }
      
      for(unsigned int i=0; i<dim; i++){
          for(unsigned int j=0; j<dim; j++)
              gsl_matrix_set(dcp.get(),i,j,0.0);
      }

    // set mixing angles
    gsl_matrix_set(th.get(),0,1,th12);
    gsl_matrix_set(th.get(),0,2,th13);
    gsl_matrix_set(th.get(),1,2,th23);
    

    gsl_matrix_complex* U = gsl_matrix_complex_alloc(dim,dim);
    gsl_matrix_complex* R = gsl_matrix_complex_alloc(dim,dim);
    gsl_matrix_complex* dummy = gsl_matrix_complex_alloc(dim,dim);
    gsl_matrix_complex_set_identity(U);
    gsl_matrix_complex_set_identity(R);
    gsl_matrix_complex_set_zero(dummy);
    
    const auto unit=gsl_complex_rect(1,0);
    const auto zero=gsl_complex_rect(0,0);
    
    //construct each subspace rotation and accumulate the product
    for(size_t j=1; j<dim; j++){
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


  // Returns the 6-dimensional transformation matrix within our theory.
  std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex *)> nuSQUIDS_SM_Copies::GetTransformationMatrix_SM_Copies() {

    // Define the PMNS matrix
    std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> PMNS = GetPMNS();
    // Calculate cos and sin
    double cos_value = std::sqrt((N - 1.0) / N);
    double sin_value = 1.0 / std::sqrt(N);

    const auto cos_gsl = gsl_complex_rect(cos_value,0);
    const auto sin_gsl = gsl_complex_rect(sin_value,0);

    // Calculate U using Kronecker product
    gsl_matrix_complex* U = gsl_matrix_complex_alloc(6, 6);
    gsl_complex element1;
    gsl_complex element2;
    for(unsigned int n=0; n<2; n++){
      for(unsigned int i=0; i<dim; i++){
        for(unsigned int j=0; j<dim; j++){
            
          element1 = gsl_complex_mul(cos_gsl, gsl_matrix_complex_get(PMNS.get(),i, j));
          element2 = gsl_complex_mul(sin_gsl, gsl_matrix_complex_get(PMNS.get(),i, j));
          element2 = gsl_complex_mul(gsl_complex_exp(gsl_complex_mul(gsl_complex_rect(n+1, 0), gsl_complex_log(gsl_complex_rect(-1, 0)))), element2);
          gsl_matrix_complex_set(U,i+3*n,j+3*n,element1);
          gsl_matrix_complex_set(U,i+3-3*n,j+3*n,element2);
        };
      };
    };


    // Wrap U in unique_ptr with custom deleter
    auto deleter = [](gsl_matrix_complex* mat) { gsl_matrix_complex_free(mat); };
    return std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex *)>(U, deleter);
  }


  void nuSQUIDS_SM_Copies::iniProjectors(){

    
    std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> U_SM_Copies = GetTransformationMatrix_SM_Copies();
    gsl_matrix_complex* U = U_SM_Copies.get();

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
        b1_proj[rho][flv] = b1_proj[rho][flv].Rotate(U);
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
          evol_b1_proj[rho][flv][e1] = evol_b1_proj[rho][flv][e1].Rotate(U);
          AntineutrinoCPFix(rho);
        }
      }
    }
  } 


  void nuSQUIDS_SM_Copies::SetIniFlavorProyectors(){

    std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> U_SM_Copies = GetTransformationMatrix_SM_Copies();
    gsl_matrix_complex* U = U_SM_Copies.get();

    for(unsigned int rho = 0; rho < nrhos; rho++){
      for(unsigned int flv = 0; flv < numneu; flv++){
        for(unsigned int e1 = 0; e1 < ne; e1++){
          evol_b1_proj[rho][flv][e1] = b0_proj[flv];

          AntineutrinoCPFix(rho);
          evol_b1_proj[rho][flv][e1] = evol_b1_proj[rho][flv][e1].Rotate(U);
          AntineutrinoCPFix(rho);
        }
        b1_proj[rho][flv] = b0_proj[flv];

        AntineutrinoCPFix(rho);
        b1_proj[rho][flv] = b1_proj[rho][flv].Rotate(U);
        AntineutrinoCPFix(rho);
      }
    }
  }
}


  