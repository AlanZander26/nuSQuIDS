#ifndef SM_COPIES_H
#define SM_COPIES_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>

namespace nusquids {

  class nuSQUIDS_SM_Copies: public nuSQUIDS {

    private:
      // Dimensions. Here always 3
      int dim=3;
      // SM_Copies parameters
      double N; // Number of SM copies.
      double mu; // Mass factor mH_i = mu*m_i. Here we assume mu_i=mu for all i.
      double m0; // Mass of the lightest neutrino state in (eV?? check this).
      // Masses of light neutrinos.
      double m1, m2, m3;
      // Masses of heavy neutrinos.
      double mH1, mH2, mH3;
      // Normal neutrino mass ordering hierarchy.
      bool NormalOrdering;


    public:

      nuSQUIDS_SM_Copies(marray<double,1> E_vector, double N, double mu, double m0, bool NormalOrdering, 
      unsigned int numneu = 3, NeutrinoType NT = both, bool iinteraction = false, 
      std::shared_ptr<CrossSectionLibrary> ncs = nullptr) : N(N), mu(mu), m0(m0), NormalOrdering(NormalOrdering),
      nuSQUIDS(E_vector, numneu, NT, iinteraction, ncs)
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
        Set_SquareMassDifference(1,Deltaq_m21); 
        Set_SquareMassDifference(2,Deltaq_m31); 

        if (NormalOrdering) {
          m1 = m0;
          m2 = std::sqrt(Deltaq_m21 + std::pow(m0, 2));
          m3 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
        } else {
          m1 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
          m2 = std::sqrt(Deltaq_m21 + Deltaq_m31 + std::pow(m0, 2));
          m3 = m0;
        }   
          
        mH1 = mu*m1;
        mH2 = mu*m2;
        mH3 = mu*m3;

        double Deltaq_mH1 = std::pow(mH1, 2) - std::pow(m0, 2);
        double Deltaq_mH2 = std::pow(mH2, 2) - std::pow(m0, 2);
        double Deltaq_mH3 = std::pow(mH3, 2) - std::pow(m0, 2);
          
        Set_SquareMassDifference(3,Deltaq_mH1); // dm^2_H1
        Set_SquareMassDifference(4,Deltaq_mH2); // dm^2_H2
        Set_SquareMassDifference(5,Deltaq_mH3); // dm^2_H3
      }

      std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetPMNS(double th12 = 0.563942, double th13 = 0.154085, double th23 = 0.785398);

      std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex *)> GetTransformationMatrix_SM_Copies();

      void iniProjectors() override;

      void SetIniFlavorProyectors() override;
      
      
      
    
    
  };
} // close nusquids namespace

#endif // SM_COPIES_H