#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include "ADD/ADD.h" // Change for another model.

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
  std::string output_path;
  if(argc != 7){
      printf("ERROR:USAGE: the amount of arguments must be 6. \n");
      exit(0);
  } else {
      input_flux_path  = argv[1];
      input_earth_path = argv[2];
      output_path      = argv[3];
      a                = atof(argv[4]);
      m0               = atof(argv[5]);
      NormalOrdering   = argv[6];

  }
  if(output_path[output_path.length()-1]!='/'){
    output_path = output_path +'/';}

  std::cout<<"Inpath Flux: "<<input_flux_path<<std::endl;
  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;

  const squids::Const units;
  nuSQUIDSAtm<nuSQUIDS_ADD> nus_atm_prompt(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350), N_KK, a, m0, NormalOrdering, neutrino_type, iinteraction);

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_earth_path); 
  nus_atm_prompt.Set_EarthModel(earth);

  nus_atm_prompt.Set_TauRegeneration(true);

  nus_atm_prompt.Set_ProgressBar(false);

  double error = 1.0e-16;
  // setup integration settings
  nus_atm_prompt.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm_prompt.Set_rel_error(error);
  nus_atm_prompt.Set_abs_error(error);

  marray<double,2> input_flux = quickread(input_flux_path);

  marray<double,4> inistate_prompt{nus_atm_prompt.GetNumCos(),nus_atm_prompt.GetNumE(),2,numneu};
  std::fill(inistate_prompt.begin(),inistate_prompt.end(),0);

  marray<double,1> cos_range = nus_atm_prompt.GetCosthRange();
  marray<double,1> e_range = nus_atm_prompt.GetERange();

  assert( input_flux.extent(0) == nus_atm_prompt.GetNumCos()*nus_atm_prompt.GetNumE() );

  for ( int ci = 0 ; ci < nus_atm_prompt.GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nus_atm_prompt.GetNumE(); ei++){
      double enu = e_range[ei]/units.GeV;
      assert( std::fabs(enu - input_flux[ci*e_range.size() + ei][1] ) < 1.e-4 );
      double cth = cos_range[ci];
      assert( std::fabs(cth- input_flux[ci*e_range.size() + ei][0] ) < 1.e-4 );

      inistate_prompt[ci][ei][0][0] = input_flux[ci*e_range.size() + ei][2];
      inistate_prompt[ci][ei][0][1] = input_flux[ci*e_range.size() + ei][4];
      inistate_prompt[ci][ei][0][2] = input_flux[ci*e_range.size() + ei][6];
      inistate_prompt[ci][ei][0][3] = 0.;

      inistate_prompt[ci][ei][1][0] = input_flux[ci*e_range.size() + ei][3];
      inistate_prompt[ci][ei][1][1] = input_flux[ci*e_range.size() + ei][5];
      inistate_prompt[ci][ei][1][2] = input_flux[ci*e_range.size() + ei][7];
      inistate_prompt[ci][ei][1][3] = 0.;
    }
  }

  nus_atm_prompt.Set_initial_state(inistate_prompt,flavor);

  nus_atm_prompt.EvolveState();

   for ( int ci = 0 ; ci < nus_atm_prompt.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm_prompt.GetNumE(); ei++){
       double enu = e_range[ei];
       double cth = cos_range[ci];
       for(unsigned int flv = 0; flv < 3; flv++){
         if(nus_atm_prompt.EvalFlavor(flv,cth,enu,0) < 0)
           std::cout << "neg nu    propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm_prompt.EvalFlavor(flv,cth,enu,0) << std::endl;
         if(nus_atm_prompt.EvalFlavor(flv,cth,enu,1) < 0)
           std::cout << "neg nubar propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm_prompt.EvalFlavor(flv,cth,enu,1) << std::endl;
       }
     }
   }

  nus_atm_prompt.WriteStateHDF5(output_path+"/prompt_atmospheric_"+"ADD_"+std::to_string(a)+"_"+std::to_string(m0)+".hdf5");
  
  std::cout << "finish" << std::endl;

  return 0;
}