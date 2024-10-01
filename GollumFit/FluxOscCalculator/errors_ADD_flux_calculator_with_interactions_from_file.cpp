#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include "ADD/ADD.h" // Change for another model.
// The following command is still not working for the compilation: gcc conventional_atmospheric_flux.cpp -I/usr/include/hdf5/serial -I../../TeVSGT/ADD -L/usr/lib/x86_64-linux-gnu -lhdf5_serial -lhdf5_cpp -lpthread

using namespace nusquids;

int main(int argc, char* argv[]){
  // parameters of the theory
  double a, m0;
  bool NormalOrdering;
  NeutrinoType neutrino_type;
  bool iinteraction = true; // Otherwise I get: 'std::runtime_error' what():  nuSQUIDS::Error::nuSQuIDs has been initialized without interactions, thus tau regeneration cannot be enabled.
  unsigned int N_KK;
  unsigned int numneu;
  std::string input_flux_dir, input_earth_path;
  std::string flux_name;
  std::string output_path;
  if(argc != 8){
      printf("ERROR:USAGE: the amount of arguments must be 7. \n");
      exit(0);
  } else {
      input_flux_dir  = argv[1];
      input_earth_path = argv[2];
      output_path      = argv[3];
      flux_name        = argv[4];
      a                = atof(argv[5]);
      m0               = atof(argv[6]);
      NormalOrdering   = argv[7];
      neutrino_type = both; // If I do not set "both", I get: 'std::runtime_error' what():  nuSQUIDS::Error::Cannot set TauRegeneration to True when NT != 'both'.//neutrino; // Change this to accept also antineutrino or both. 
      N_KK             = 2;
      numneu = 3*(N_KK+1);

  }
  if(output_path[output_path.length()-1]!='/'){
    output_path = output_path +'/';}

  std::cout<<"Indir Flux: "<<input_flux_dir<<std::endl;
  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;
  std::cout<<"Flux name: "<<flux_name<<std::endl;
    
  const squids::Const units;
  nuSQUIDSAtm<nuSQUIDS_ADD> nus_atm(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350), N_KK, a, m0, NormalOrdering, neutrino_type, iinteraction);

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
  marray<double,2> input_flux = quickread(input_flux_dir+"/ddm_"+flux_name+".dat"); 

  // construct the kaon initial state
   marray<double,4> inistate {nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
   std::fill(inistate.begin(),inistate.end(),0);

   marray<double,1> cos_range = nus_atm.GetCosthRange();
   marray<double,1> e_range = nus_atm.GetERange();

   assert( input_flux.extent(0) == nus_atm.GetNumCos()*nus_atm.GetNumE() );

   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei]/units.GeV;
       assert( std::fabs(enu - input_flux[ci*e_range.size() + ei][1] ) < 1.e-4 );
       double cth = cos_range[ci];
       assert( std::fabs(cth- input_flux[ci*e_range.size() + ei][0] ) < 1.e-4 );

       inistate[ci][ei][0][0] = input_flux[ci*e_range.size() + ei][2];
       inistate[ci][ei][0][1] = input_flux[ci*e_range.size() + ei][4];
       inistate[ci][ei][0][2] = input_flux[ci*e_range.size() + ei][6];
       inistate[ci][ei][0][3] = 0.;

       inistate[ci][ei][1][0] = input_flux[ci*e_range.size() + ei][3];
       inistate[ci][ei][1][1] = input_flux[ci*e_range.size() + ei][5];
       inistate[ci][ei][1][2] = input_flux[ci*e_range.size() + ei][7];
       inistate[ci][ei][1][3] = 0.;
     }
   }

  nus_atm.Set_initial_state(inistate,flavor);

  nus_atm.EvolveState();

   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei];
       double cth = cos_range[ci];
       for(unsigned int flv = 0; flv < 3; flv++){
         if(nus_atm.EvalFlavor(flv,cth,enu,0) < 0)
           std::cout << "neg nu    propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm.EvalFlavor(flv,cth,enu,0) << std::endl;
         if(nus_atm.EvalFlavor(flv,cth,enu,1) < 0)
           std::cout << "neg nubar propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm.EvalFlavor(flv,cth,enu,1) << std::endl;
       }
     }
   }

  nus_atm.WriteStateHDF5(output_path+"/"+flux_name+"_ADD_"+std::to_string(a)+"_"+std::to_string(m0)+".hdf5");

  std::cout << "finish" << std::endl;

  return 0;
}
