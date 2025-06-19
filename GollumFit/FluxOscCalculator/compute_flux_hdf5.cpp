#include <vector>
#include <iostream>
#include <fstream>

#include <nuSQuIDS/nuSQuIDS.h>

#include "ADD/ADD.h"
#include "SM_Copies/SM_Copies.h"

using namespace nusquids;

int main(int argc, char* argv[]){

    std::string nusquids_data_file_path;
    nusquids_data_file_path = argv[1];
    
    // Find the last path separator (either / or \ for cross-platform support)
    size_t pos = nusquids_data_file_path.find_last_of("/\\");

    // Extract directory path
    std::string path_to_dir = (pos != std::string::npos) ?
                            nusquids_data_file_path.substr(0, pos) :
                            "";

    bool atmospheric_height_randomization = false;

    const squids::Const units;

    double E_min = 1.e2;
    double E_max = 1.e6;
    double czmin=-1.;
    double czmax=0.2;

    nuSQUIDSAtm<nuSQUIDS_SM_Copies> nus_add_atm(nusquids_data_file_path);
    nuSQUIDSAtm<> nus_atm(nusquids_data_file_path);

    unsigned int numneu = nus_add_atm.GetNumNeu();  

    std::ofstream file_d(path_to_dir+"/"+"flux_final_hdf5_derived.txt");

    int Nen=350;
    int Ncz=100;
    double lEmin=log10(E_min*units.GeV);
    double lEmax=log10(E_max*units.GeV);
    double dcz = (czmax - czmin) / (Ncz - 1);
    double dLE = (lEmax - lEmin) / (Nen - 1);

    //Writing to the file!  
    file_d << "# log10(E) cos(zenith) flux_i . . . ." << std::endl;
    for(int ei = 0; ei < Nen; ++ei) {
        double lE = lEmin + ei * dLE;
        double E = pow(10.0, lE);
    
        for(int czi = 0; czi < Ncz; ++czi) {
        double cz = czmin + czi * dcz;
    
        file_d << lE << " " << cz;
        for(int fl = 0; fl < numneu; ++fl) {
            for(int rho = 0; rho < 2; ++rho) {
            file_d << " " << nus_add_atm.EvalFlavor(fl, cz, E, rho);
            }
        }
        file_d << std::endl;
        }
        file_d << std::endl;
    }


    std::ofstream file_b(path_to_dir+"/"+"flux_final_hdf5_base.txt");

    //Writing to the file!  
    file_b << "# log10(E) cos(zenith) flux_i . . . ." << std::endl;
    for(int ei = 0; ei < Nen; ++ei) {
        double lE = lEmin + ei * dLE;
        double E = pow(10.0, lE);
    
        for(int czi = 0; czi < Ncz; ++czi) {
        double cz = czmin + czi * dcz;
    
        file_b << lE << " " << cz;
        for(int fl = 0; fl < numneu; ++fl) {
            for(int rho = 0; rho < 2; ++rho) {
            file_b << " " << nus_atm.EvalFlavor(fl, cz, E, rho);
            }
        }
        file_b << std::endl;
        }
        file_b << std::endl;
    }



    return 0;
}