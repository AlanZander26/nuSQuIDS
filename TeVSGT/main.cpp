#define USE_SM_COPIES
// USE_ADD (ADD), USE_SM_COPIES (SM_Copies), USE_SM (SM)
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <nuSQuIDS/nuSQuIDS.h>

#ifdef USE_ADD
#include "ADD/ADD.h" // Modify this to absolute paths

#elif defined(USE_SM_COPIES)
#include "SM_Copies/SM_Copies.h"

#elif defined(USE_SM)
#else 
#error "Error: Model type not defined."

#endif


using namespace nusquids;

int main () {

  squids::Const units;

#ifdef USE_ADD
  double a; // Radius of largest extra dimension [micro m]
  double m0; // Mass of lightest neutrino in the SM [eV]
  std::cin >> a >> m0;
  unsigned int N_KK = 2;
  unsigned int numneu = 3*(N_KK + 1); 

#elif defined(USE_SM_COPIES)
  double N; // Number of SM copies
  double mu; // Mass factor
  double m0; // Mass of lightest neutrino in the SM [eV]
  std::cin >> N >> mu >> m0;
  unsigned int numneu = 6;

#elif defined(USE_SM)
  unsigned int numneu = 3;
#endif

  double E_min, E_max; // Min. and max. energy, respectively [in GeV]
  std::string medium; // Medium of propagation ["vacuum", "earth"]
  double medium_param_min;
  double medium_param_max;
  unsigned int N_medium_param;
  std::string neutrino_type_str;
  NeutrinoType neutrino_type; // Neutrino type ["neutrino", "antineutrino"]
  std::string NormalOrdering_str;
  bool NormalOrdering; // Ordering of the masses
  int N_energy_grid = 100; // Number of the internal energy grid points of the nuSQuIDS object
  int Nen; // Number of energies we want the result
  std::cin >> E_min >> E_max >> medium >> medium_param_min >> medium_param_max >> N_medium_param
  >> neutrino_type_str >> NormalOrdering_str >> N_energy_grid >> Nen;

  if (neutrino_type_str == "neutrino") {
    neutrino_type = neutrino;
  }
  else if (neutrino_type_str == "antineutrino") {
    neutrino_type = antineutrino;
  }
  else {
    std::cout << "Error: neutrino_type must be either 'neutrino' or 'antineutrino'." << std::endl;
    return 1;
  }

  if (NormalOrdering_str == "true") {
    NormalOrdering = true;
  }
  else if (NormalOrdering_str == "false") {
    NormalOrdering = false;
  }
  else {
    std::cout << "Error: NormalOrdering must be either 'true' or 'false'." << std::endl;
    return 1;
  }

  double initial_flux_ratio; // Initial flavor flux ratios. The components need to add up to 1.
  std::vector<double> initial_flux_ratios(3);
  double sum_ratios=0.;
  for(int i = 0 ; i < 3; i++){
    std::cin >> initial_flux_ratio;
    initial_flux_ratios[i] = initial_flux_ratio;
    sum_ratios+=initial_flux_ratio;
  }
  
  if(sum_ratios != 1.){
    std::cout << "Error: the initial flavor flux ratios must add to 1." << std::endl;
    return 1;
  }


  std::vector<double> medium_param_vec = linspace_vec(medium_param_min, medium_param_max, N_medium_param);
  double medium_param;

  for(int i = 0; i < N_medium_param; ++i) {
    medium_param = medium_param_vec[i];

#ifdef USE_ADD
    nuSQUIDS_ADD nus(logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), N_KK, a, m0, NormalOrdering, neutrino_type);

#elif defined(USE_SM_COPIES)
    nuSQUIDS_SM_Copies nus(logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), N, mu, m0, NormalOrdering, numneu, neutrino_type);

#elif defined(USE_SM)
    nuSQUIDS nus(logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), numneu, neutrino_type, false);

#endif

    if (medium == "vacuum") {
      double L; // Baseline [km]
      L = medium_param;

      const float L_baseline = L*units.km;
      std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
      std::shared_ptr<Vacuum::Track> track_vacuum = std::make_shared<Vacuum::Track>(L_baseline);
      nus.Set_Body(vacuum);
      nus.Set_Track(track_vacuum);
  
  } else if (medium == "earth") {
      double phi; // Zenit angle. phi = acos(-1.) when neutrinos crossing the Earth.
      phi = std::acos(medium_param);
      //Here we define the trajectory that the particle follows and the object for more examples
      // of how construct a track and object look body_track example.
      //Declaration of the body, EarthAtm is one of the predefined bodies
      std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
      //Definition of the track, in encodes the trajectory inside the body, here is declared with the zenith angle.
      auto track_atm = std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(phi));
      //We set this in the nusSQuID object.
      nus.Set_Body(earth_atm);
      nus.Set_Track(track_atm);
  } else {
      std::cout << "Error: invalid medium of propagation." << std::endl;
  }
    
    //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
    nus.Set_h_max( 500.0*units.km );
    nus.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus.Set_rel_error(1.0e-5);
    nus.Set_abs_error(1.0e-5);

    //Construct the initial state

    //E_range is an array that contains all the energies.
    marray<double,1> E_range = nus.GetERange();
    //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
    marray<double,2> inistate{E_range.size(),numneu};
    double N0 = 1.0e18;
    // Set a power-law spectra for the neutrinos of the SM with initial flavor flux ratios given in "initial_flux_ratios". BSM iniital flux flavors are set to zero. To change this, change the code below.
    for ( int i = 0 ; i < inistate.extent(0); i++){
        for ( int k = 0; k < 3; k ++){
          inistate[i][k] = initial_flux_ratios[k]*N0*pow(E_range[i],-2);
        }
        for ( int k = 3; k < inistate.extent(1); k ++){
          inistate[i][k] = 0.0;
        }
    } 
    nus.Set_initial_state(inistate,flavor);
    nus.EvolveState();

    // Creating a dynamic file name based on parameters
    std::ostringstream filename;
    filename << std::fixed << std::setprecision(2);
#ifdef USE_ADD
    filename << "output_ADD_a_" << a << "_m0_" << m0 << "_";
#elif defined(USE_SM_COPIES)
    filename << "output_SM_Copies_N_" << N << "_mu_" << mu << "_m0_" << m0 << "_";
#elif defined(USE_SM)
    filename << "output_SM_";
#endif
    filename << "Emin_" << E_min << "_Emax_" << E_max << "_" << medium << "_";
    filename << "Nengrid_" << N_energy_grid << ".hdf5";

    nus.WriteStateHDF5(filename.str());
  }
  std::cout << std::endl <<  "Done! Merging HDF5 files..." << std::endl;    
  const char* command = "python3 combine_hdf5.py";
  int result = std::system(command);
  if (result == 0) {
        std::cout << "Done!." << std::endl;
  } else {
        std::cerr << "Failed to merge HDF5 files." << std::endl;
  }
  return 0;
}