#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <LeptonWeighter/Flux.h>
#include <LeptonWeighter/nuSQFluxInterface.h>
#include "ADD/ADD.h"

using namespace nusquids;

int main(int argc, char* argv[]){
  // parameters of the theory
  double a, m0;
  bool NormalOrdering;
  NeutrinoType neutrino_type;
  bool iinteraction = true; // Otherwise we get: 'std::runtime_error' what():  nuSQUIDS::Error::nuSQuIDs has been initialized without interactions, thus tau regeneration cannot be enabled.
  neutrino_type = both; // If we do not set "both", we get: 'std::runtime_error' what():  nuSQUIDS::Error::Cannot set TauRegeneration to True when NT != 'both'.//neutrino; // Change this to accept also antineutrino or both. 
  unsigned int N_KK = 2;
  unsigned int numneu = 3*(N_KK+1);
  std::string input_flux_path, input_earth_path;
  std::string output_path, flux_name;
  if(argc != 8){
      printf("ERROR:USAGE: the amount of arguments must be 7. \n");
      exit(0);
  } else {
      input_flux_path  = argv[1]; // This path is not needed for this script but is kept to mantain the same format as the other cpp files.
      input_earth_path = argv[2];
      output_path      = argv[3];
      flux_name        = argv[4]; // The flux name is not needed for this script but is kept to mantain the same format as the other cpp files. >
      a                = atof(argv[5]);
      m0               = atof(argv[6]);
      NormalOrdering   = argv[7];

  }
  if(output_path[output_path.length()-1]!='/'){
    output_path = output_path +'/';}

  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;

  double baseline_astro_normalization = 1.0e-18; // nu/GeV/s/cm^2/sr
  double baseline_astro_spectral_index = -2.5; // center of things
  auto fluxAstro_ = std::make_shared<LW::PowerLawFlux>(baseline_astro_normalization,
                                                       baseline_astro_spectral_index);

  const squids::Const units;
  nuSQUIDSAtm<nuSQUIDS_ADD> nus_atm_astro(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350), N_KK, a, m0, NormalOrdering, neutrino_type, iinteraction);

  nus_atm_astro.Set_TauRegeneration(true);

  nus_atm_astro.Set_IncludeOscillations(true);
  nus_atm_astro.Set_ProgressBar(false);

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_earth_path);
  nus_atm_astro.Set_EarthModel(earth);

  double error = 1.0e-17;
  // setup integration settings
  nus_atm_astro.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm_astro.Set_rel_error(error);
  nus_atm_astro.Set_abs_error(error);


  // construct the astro initial state
   marray<double,4> inistate_astro {nus_atm_astro.GetNumCos(),nus_atm_astro.GetNumE(),2,numneu};
   std::fill(inistate_astro.begin(),inistate_astro.end(),0);

   marray<double,1> cos_range = nus_atm_astro.GetCosthRange();
   marray<double,1> e_range = nus_atm_astro.GetERange();

  LW::Event scratch_lw_e;
   for ( int ci = 0 ; ci < nus_atm_astro.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm_astro.GetNumE(); ei++){
       double enu = e_range[ei]/units.GeV;
       double cth = cos_range[ci];

       scratch_lw_e.energy=enu;
       scratch_lw_e.zenith=acos(cth);

       inistate_astro[ci][ei][0][0] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][0][1] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][0][2] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][0][3] = 0.;

       inistate_astro[ci][ei][1][0] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][1][1] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][1][2] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][1][3] = 0.;
     }
   }

  nus_atm_astro.Set_initial_state(inistate_astro,flavor);

  nus_atm_astro.EvolveState();


   for ( int ci = 0 ; ci < nus_atm_astro.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm_astro.GetNumE(); ei++){
       double enu = e_range[ei];
       double cth = cos_range[ci];
       for(unsigned int flv = 0; flv < 3; flv++){
         if(nus_atm_astro.EvalFlavor(flv,cth,enu,0) < 0)
           std::cout << "neg nu    propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm_astro.EvalFlavor(flv,cth,enu,0) << std::endl;
         if(nus_atm_astro.EvalFlavor(flv,cth,enu,1) < 0)
           std::cout << "neg nubar propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm_astro.EvalFlavor(flv,cth,enu,1) << std::endl;
       }
     }
   }


  nus_atm_astro.WriteStateHDF5(output_path+"/astro_"+"ADD_"+std::to_string(a)+"_"+std::to_string(m0)+".hdf5");

  std::cout << "finish" << std::endl;

  return 0;
}
