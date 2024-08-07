#define USE_SM
// USE_ADD (ADD), USE_SM_COPIES (SM_Copies), USE_SM (SM)
#include <vector>
#include <iostream>
#include <fstream>
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


std::vector<double> linspace_vec(double a,double b,unsigned int N){
  if(N==0)
      throw std::length_error("number of samples requested from linspace must be nonzero");
  std::vector<double> linpoints(N);
  if(N==1){
      if(a==b){
          linpoints[0] = a;
          return linpoints;
      }
      else{
          throw std::invalid_argument("When N = 1, a and b must be equal.");
      }
  }
  double step_lin = (b - a)/double(N-1);
  
  double c = a;
  for(unsigned int i=0; i<N-1; i++, c+=step_lin)
      linpoints[i] = c;
  linpoints[N-1] = b;

  return linpoints;
}


int main () {

  squids::Const units;

#ifdef USE_ADD
  double a_min, a_max; // Radius range of largest extra dimension [micro m]
  unsigned int N_a; // # of bins in dimension "a"
  double m0_min, m0_max; // Mass range of lightest neutrino in the SM [eV]
  unsigned int N_m0; // # of bins in dimension "m0"
  std::cin >> a_min >> a_max >> N_a >> m0_min >> m0_max >> N_m0;
  unsigned int N_KK = 2;
  unsigned int numneu = 3*(N_KK + 1); 

  std::vector<double> a_vec = linspace_vec(a_min, a_max, N_a);
  std::vector<double> m0_vec = linspace_vec(m0_min, m0_max, N_m0);
  double a;
  double m0;
#elif defined(USE_SM_COPIES)
  double N_min, N_max; // Range of number of SM copies
  unsigned int N_N; // # of bins in dimension "N"
  double mu_min, mu_max; // Range of the mass factor
  unsigned int N_mu; // # of bins in dimension "mu"
  double m0_min, m0_max; // Mass range of lightest neutrino in the SM [eV]
  unsigned int N_m0; // # of bins in dimension "m0"
  std::cin >> N_min >> N_max >> N_N >> mu_min >> mu_max >> N_mu >> m0_min >> m0_max >> N_m0;
  unsigned int numneu = 6;

  std::vector<double> N_vec = linspace_vec(N_min, N_max, N_N);
  std::vector<double> mu_vec = linspace_vec(mu_min, mu_max, N_mu);
  std::vector<double> m0_vec = linspace_vec(m0_min, m0_max, N_m0);
  double N;
  double mu;
  double m0;
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
  int N_energy_grid = 200; // Number of the internal energy grid points of the nuSQuIDS object
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
  //In this part we will save the values in a txt file to be able to plot or manipulate later.
  //Notice that this is not going to have all the information about the quantum evolution, for that 
  //we need to save the information using the HDF5 writing function.
  std::ofstream file("fluxes_flavor.txt");

  for(int i = 0; i < N_medium_param; ++i) {
    medium_param = medium_param_vec[i];

#ifdef USE_ADD
    for(int j = 0; j < N_m0; ++j) {
        m0 = m0_vec[j];
        for(int k = 0; k < N_a; ++k) {
            a = a_vec[k];
            nuSQUIDS_ADD nus(logspace(E_min*units.GeV,E_max*units.GeV,N_energy_grid), N_KK, a, m0, NormalOrdering, neutrino_type);

#elif defined(USE_SM_COPIES)
    for(int j = 0; j < N_m0; ++j) {
        m0 = m0_vec[j];
        for(int k = 0; k < N_mu; ++k) {
            mu = mu_vec[k];
            for(int l = 0; l < N_N; ++l) {
                N = N_vec[l];
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

          //We set the GSL step function
          nus.Set_GSL_step(gsl_odeiv2_step_rk4);

          //Setting the numerical precision of gsl integrator.
          nus.Set_rel_error(1.0e-5);
          nus.Set_abs_error(1.0e-5);

          //Set true the progress bar during the evolution.
          //nus.Set_ProgressBar(true);

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


          //Set the initial state in nuSQuIDS object
          nus.Set_initial_state(inistate,flavor);
          //Propagate the neutrinos in the earth for the path defined in path
          nus.EvolveState();
        
          //number of energies we want the result, notice that this can be larger than the number of the internal grid of 
          //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
          //and vacuum oscillations are solved analytically for the given energy.
          double lEmin=std::log10(E_min);
          double lEmax=std::log10(E_max);
          std::vector<double> lE_vec = linspace_vec(lEmin, lEmax, Nen);
          double lE;

          #ifdef USE_ADD
          file << "# Energy[GeV] Prob_e Prob_mu Prob_tau " << "(a [um] = " << a << ", m0 [eV] = " << m0 << ", medium_param = " << medium_param << ")" << std::endl;

          #elif defined(USE_SM_COPIES)
          file << "# Energy[GeV] Prob_e Prob_mu Prob_tau " << "(N = " << N << ", mu = " << mu << ", m0 [eV] = " << m0 << ", medium_param = " << medium_param << ")" << std::endl;

          #elif defined(USE_SM)
          file << "# Energy[GeV] Prob_e Prob_mu Prob_tau " << "(medium_param = " << medium_param << ")" << std::endl;

          #endif
          for(double i=0; i<Nen; i++){
            lE = lE_vec[i];
            double E=pow(10.0,lE); // This expression is dimensionless of the form energy/GeV
            file << E << " ";
            for(int fl=0; fl<3; fl++){
              file << " " <<  nus.EvalFlavor(fl, E*units.GeV)/(N0*pow(E*units.GeV,-2));
            }
            file << std::endl;
          }
        
#ifdef USE_ADD
        }
    }
#elif defined(USE_SM_COPIES)
      }
    }
  }
#elif defined(USE_SM)
#endif
  }
  file.close();
  std::cout << std::endl <<  "Done! " << std::endl;    
  

  return 0;

}



