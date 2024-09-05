#ifndef DarkDim_H
#define DarkDim_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids {

    class nuSQUIDS_DarkDim: public nuSQUIDS {

        private:
            unsigned int N_KK; // Number of KK states
            unsigned int dim_DarkDim; // Dimension of the DarkDim space = 3*(N_KK + 1)
            double a; // Radius of extra dimension [mu m]
            double m0; // Mass of the lightest neutrino state [eV]
            std::array<double, 3> ca; // List of bare masses c_i multiplied by the radius of dark dimension [dimensionless].
            //std::array<double, 3> lambda; // Effective Yukawa couplings.
            //const double Mf; // Fundamental scale of gravity [GeV].
            double m1, m2, m3; // Masses of light neutrinos.
            
            bool NormalOrdering; // Normal neutrino mass ordering hierarchy.

            gsl_vector *Lambdaq; // Vector with "squared masses" (lambdas, see paper Machado et.al.) in ascending order.
            gsl_matrix_complex *W; // Transformation matrix between flavor and mass basis.

        public:

            nuSQUIDS_DarkDim(marray<double,1> E_vector, unsigned int N_KK, double a, double m0, std::array<double, 3> ca, 
            bool NormalOrdering, NeutrinoType NT = both, bool iinteraction = false, 
            std::shared_ptr<CrossSectionLibrary> ncs = nullptr) : N_KK(N_KK), dim_DarkDim(3*(N_KK + 1)), a(a), m0(m0), ca(ca), NormalOrdering(NormalOrdering), 
            nuSQUIDS(E_vector, 3*(N_KK + 1), NT, iinteraction, ncs)
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

                Lambdaq = gsl_vector_alloc(dim_DarkDim-1);
                W = gsl_matrix_complex_alloc(dim_DarkDim, dim_DarkDim);
                iniMatrices(Lambdaq, W, th01, th02, th12);

                for (int j = 0; j < dim_DarkDim-1; j++) {
                    Set_SquareMassDifference(j+1,gsl_vector_get(Lambdaq, j));
                }

                //===============================
                // set the lambdas to reproduce the right light neutrino masses   //
                //===============================  

                //lambda = calculate_lambda();

            }

            std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetPMNS(double th12 = 0.563942, double th13 = 0.154085, double th23 = 0.785398);

            void iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th12, double th13, double th23);

            void iniProjectors() override;

            void SetIniFlavorProyectors() override;     

            //std::array<double, 3> calculate_lambda(); 
        
        
    };
} // close nusquids namespace

#endif // DarkDim_H