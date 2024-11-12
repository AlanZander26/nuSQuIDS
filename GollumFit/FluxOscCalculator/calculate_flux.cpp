#define USE_ADD
// USE_ADD (ADD), USE_SM_COPIES (SM_Copies), USE_DARKDIM (DarkDim), USE_SM (SM)
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <LeptonWeighter/Flux.h>
#include <LeptonWeighter/nuSQFluxInterface.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/nuSQuIDS.h>

#ifdef USE_ADD
#include "ADD/ADD.h"

#elif defined(USE_SM_COPIES)
#include "SM_Copies/SM_Copies.h"

#elif defined(USE_DARKDIM)
#include "DarkDim/DarkDim.h"

#elif defined(USE_SM)
#else 
#error "Error: Model type not defined."

#endif


using namespace nusquids;

int main(int argc, char* argv[]){
  unsigned int numneu;
  bool NormalOrdering;
  NeutrinoType neutrino_type;
  bool iinteraction = true; // Otherwise I get: 'std::runtime_error' what():  nuSQUIDS::Error::nuSQuIDs has been initialized without interactions, thus tau regeneration cannot be enabled.
  std::string input_flux_path, input_earth_path;
  std::string output_path, flux_type;
  flux_type        = argv[1];
  output_path      = argv[2];
  NormalOrdering   = argv[3];
  neutrino_type = both; // If I do not set "both", I get: 'std::runtime_error' what():  nuSQUIDS::Error::Cannot set TauRegeneration to True when NT != 'both'.//neutrino; // Change this to accept also antineutrino or both. 
  // parameters of the theory

  // Check if flux_type is valid
  if (flux_type != "conventional" && flux_type != "prompt" && flux_type != "astro") {
      std::cerr << "Error: Invalid flux_type \"" << flux_type << "\".\n"
                << "Allowed values are: \"conventional\", \"prompt\", or \"astro\".\n";
      return 1;
  }

  if (flux_type == "conventional" | flux_type == "prompt") {
    input_flux_path  = "Data/v0.6.0_nodeis/ddm_"+flux_type+"_bestfit.dat";
  }
  else{
    input_flux_path = "";
  }

  input_earth_path = "Data/EARTH_MODEL_PREM.dat";

  #ifdef USE_ADD
  double a, m0;
  unsigned int N_KK = 2;
  numneu = 3*(N_KK+1);
  if(argc != 6){
      printf("ERROR:USAGE: the amount of arguments for ADD must be 4. \n");
      exit(0);
  } else {
      a                = atof(argv[4]);
      m0               = atof(argv[5]);
  }

  #elif defined(USE_SM_COPIES)
  double N, mu, m0;
  numneu = 6;
  if(argc != 7){
      printf("ERROR:USAGE: the amount of arguments for SM_Copies must be 5. \n");
      exit(0);
  } else {
      N                = atof(argv[4]);
      mu               = atof(argv[5]);
      m0               = atof(argv[6]);
  }

  #elif defined(USE_SM)
  numneu = 3;
  if(argc != 4){
      printf("ERROR:USAGE: the amount of arguments for the SM must be 2. \n");
      exit(0);
  }
  #endif

  if(output_path[output_path.length()-1]!='/'){
    output_path = output_path +'/';}

  std::cout<<"Inpath Flux: "<<input_flux_path<<std::endl;
  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;

  if (flux_type == "astro") {
        double baseline_astro_normalization = 1.0e-18; // nu/GeV/s/cm^2/sr
        double baseline_astro_spectral_index = -2.5;   // center of things
        auto fluxAstro_ = std::make_shared<LW::PowerLawFlux>(
            baseline_astro_normalization,
            baseline_astro_spectral_index
        );
  }
    
  const squids::Const units;

  double E_min = 1.e2;
  double E_max = 1.e6;
  double czmin=-1.;
  double czmax=0.2;
  unsigned int N_energy_grid = 350;
  unsigned int N_cz_grid = 100;
  #ifdef USE_ADD
  nuSQUIDSAtm<nuSQUIDS_ADD> nus_atm(linspace(czmin,czmax,N_cz_grid), logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), N_KK, a, m0, NormalOrdering, neutrino_type, iinteraction);

  #elif defined(USE_SM_COPIES)
  nuSQUIDSAtm<nuSQUIDS_SM_Copies> nus_atm(linspace(czmin,czmax,N_cz_grid), logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), N, mu, m0, NormalOrdering, numneu, neutrino_type, iinteraction);

  #elif defined(USE_SM)
    nuSQUIDSAtm<> nus_atm(linspace(czmin,czmax,N_cz_grid), logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), numneu, neutrino_type, iinteraction);

  #endif

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_earth_path); 
  nus_atm.Set_EarthModel(earth);

  nus_atm.Set_TauRegeneration(true);

  nus_atm.Set_ProgressBar(false);

  double error = 1.0e-15;
  // setup integration settings
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm.Set_rel_error(error);
  nus_atm.Set_abs_error(error);


  // loading kaon and pion flux files
  marray<double,2> input_flux = quickread(input_flux_path);

  // construct the kaon initial state
   marray<double,4> inistate {nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
   std::fill(inistate.begin(),inistate.end(),0.0);

   marray<double,1> cos_range = nus_atm.GetCosthRange();
   marray<double,1> e_range = nus_atm.GetERange();

if (flux_type == "astro") {
    LW::Event scratch_lw_e;
   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei]/units.GeV;
       double cth = cos_range[ci];

       scratch_lw_e.energy=enu;
       scratch_lw_e.zenith=acos(cth);

       inistate[ci][ei][0][0] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][0][1] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][0][2] = (*fluxAstro_)(scratch_lw_e);

       inistate[ci][ei][1][0] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][1][1] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][1][2] = (*fluxAstro_)(scratch_lw_e);
     }
   }  
}
else {
  assert( input_flux.extent(0) == nus_atm.GetNumCos()*nus_atm.GetNumE() );

  // Populate only the non-zero entries in the loops
  for (int ci = 0; ci < nus_atm.GetNumCos(); ci++) {
    for (int ei = 0; ei < nus_atm.GetNumE(); ei++) {
        double enu = e_range[ei] / units.GeV;
        assert(std::fabs(enu - input_flux[ci * e_range.size() + ei][1]) < 1.e-4);
        double cth = cos_range[ci];
        assert(std::fabs(cth - input_flux[ci * e_range.size() + ei][0]) < 1.e-4);

        // Set only the relevant non-zero entries
        inistate[ci][ei][0][0] = input_flux[ci * e_range.size() + ei][2];
        inistate[ci][ei][0][1] = input_flux[ci * e_range.size() + ei][4];
        inistate[ci][ei][0][2] = input_flux[ci * e_range.size() + ei][6];
        
        inistate[ci][ei][1][0] = input_flux[ci * e_range.size() + ei][3];
        inistate[ci][ei][1][1] = input_flux[ci * e_range.size() + ei][5];
        inistate[ci][ei][1][2] = input_flux[ci * e_range.size() + ei][7];
    }
  }
}


  nus_atm.Set_initial_state(inistate,flavor);

  std::ofstream file_i(output_path+"/"+flux_type+"_flux_initial.txt");
  
  int Nen =700;
  int Ncz=100;
  double lEmin=log10(E_min*units.GeV);
  double lEmax=log10(E_max*units.GeV);

  //Writing to the file_i!  
  file_i << "# log10(E) cos(zenith) flux_i . . . ." << std::endl;
  for(double cz=czmin;cz<czmax;cz+=(czmax-czmin)/(double)Ncz){
    for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
      double E=pow(10.0,lE);
      file_i << lE << " " << cz;
      for(int fl=0; fl<numneu; fl++){
        for(int rho=0; rho<2; rho++){
	file_i << " " <<  nus_atm.EvalFlavor(fl,cz, E, rho);
      }}
      file_i << std::endl;
    }
    file_i << std::endl;
  }
  

  nus_atm.EvolveState();

   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei];
       double cth = cos_range[ci];
       for(unsigned int flv = 0; flv < numneu; flv++){
         if(nus_atm.EvalFlavor(flv,cth,enu,0) < 0)
           std::cout << "neg nu    propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm.EvalFlavor(flv,cth,enu,0) << std::endl;
         if(nus_atm.EvalFlavor(flv,cth,enu,1) < 0)
           std::cout << "neg nubar propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm.EvalFlavor(flv,cth,enu,1) << std::endl;
       }
     }
   }

nus_atm.WriteStateHDF5(output_path+"/"+flux_type+"_" + 
                       #ifdef USE_ADD
                       "ADD_" + std::to_string(a) + "_" + std::to_string(m0) +
                       #elif defined(USE_SM_COPIES)
                       "SM_Copies_" + std::to_string(N) + "_" + std::to_string(mu) + "_" + std::to_string(m0) +
                       #elif defined(USE_SM)
                       "SM" +
                       #endif
                       ".hdf5");
                       
  std::ofstream file(output_path+"/"+flux_type+"_flux_final.txt");

  //int Nen =700;
  //int Ncz=100;
  //double lEmin=log10(E_min*units.GeV);
  //double lEmax=log10(E_max*units.GeV);

  //Writing to the file!  
  file << "# log10(E) cos(zenith) flux_i . . . ." << std::endl;
  for(double cz=czmin;cz<czmax;cz+=(czmax-czmin)/(double)Ncz){
    for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
      double E=pow(10.0,lE);
      file << lE << " " << cz;
      for(int fl=0; fl<numneu; fl++){
        for(int rho=0; rho<2; rho++){
	file << " " <<  nus_atm.EvalFlavor(fl,cz, E, rho);
      }}
      file << std::endl;
    }
    file << std::endl;
  }

  std::cout << "finish" << std::endl;

  return 0;
}
