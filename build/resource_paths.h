#include <cstring>
#define SOURCE_DATA_PATH "/home/agurruty/Projects/MasterThesis/nuSQuIDS/data/"
#define INSTALL_DATA_PATH "/usr/local/share/nuSQuIDS/"

namespace{

const char dataPathBitOffValue='A';
//This needs to be volatile so that the compiler can't optimize away reads from it,
//which is important since we plan to change the value of the constant between compile-time 
//and run-time
const volatile char dataPathBitContainer[] = "NUSQUIDS_DATA_PATH_INSTALL_FLAG_PREFIX_A_NUSQUIDS_DATA_PATH_INSTALL_FLAG_SUFFIX";

bool getInstallBit(){
	static const std::size_t prefixLen = strlen("NUSQUIDS_DATA_PATH_INSTALL_FLAG_PREFIX_");
	return dataPathBitContainer[prefixLen]!='A';
}

}


