#include <string>

#ifndef ADDRESS_H
#define ADDRESS_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <iomanip>

namespace main_ns
{
#include <string>
namespace address_ns
{
#include <string>

class address_cls {
#include <string>


  private:
    string TempS;         // temporary variable for reading strings from input files

	public:
    // variables
    string Name;          // name of the input file
    string Directory;     // Input/output directory

    string Input_Dir;            // Input directory
    string OutputMatlab_Dir;     // directory to write the input file for Matlab visualizer interface
    string Info_Dir;             // directory to write general information
    string FullFile_Dir;         // directory to write the full results in the time domain analysis
    string HistoryFile_Dir;      // directory to write the time history of displacement in the time domain analysis
    string TransferFunction_Dir; // directory to write the frequency domain results

    // functions
  	address_cls();
		~address_cls();
    void address_fn();
	
};
}
}
#endif
