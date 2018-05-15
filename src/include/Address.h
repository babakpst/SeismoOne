

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <iomanip>

#ifndef ADDRESS_H
#define ADDRESS_H

namespace main_ns
{

namespace address_ns
{


class address_cls {

  private:
    std::string TempS;         // temporary variable for reading strings from input files

	public:
    // variables
    std::string Name;          // name of the input file
    std::string Directory;     // Input/output directory

    std::string Input_Dir;            // Input directory
    std::string OutputMatlab_Dir;     // directory to write the input file for Matlab visualizer interface
    std::string Info_Dir;             // directory to write general information
    std::string FullFile_Dir;         // directory to write the full results in the time domain analysis
    std::string HistoryFile_Dir;      // directory to write the time history of displacement in the time domain analysis
    std::string TransferFunction_Dir; // directory to write the frequency domain results

    // functions
  	address_cls();
    void address_fn();
	
};
}
}
#endif
