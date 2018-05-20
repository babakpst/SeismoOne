
#include "../include/Address.h"

// Constructor
main_ns::address_ns::address_cls::address_cls(){
std::cout << " --------------- SeismoOne ---------------" << std::endl;
std::cout << " ----- Developed by: Babak Poursartip ----" << std::endl;
std::cout << " -------------- version 3.0 --------------" << std::endl;
std::cout << std::endl;
}

/*
##################################################################################################

Purpose: This function reads the data file name and directories

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 05/14/2018 - Initiated: Compiled without error for the first time.

##################################################################################################
*/
void main_ns::address_ns::address_cls::address_fn(){

// Open address files 
std::cout << " opening the address file ..." << std::endl;
std::ifstream Addressfile;
Addressfile.open ("Address.txt", std::ios::in );

// Reading simulation information
std::cout << " Reading file information ..." << std::endl;
getline (Addressfile,TempS);
Addressfile >> Name ; // Input file name

getline (Addressfile,TempS);
getline (Addressfile,TempS);
Addressfile >> Directory ; // Directory of the input file

std::cout << std::endl;
std::cout << " Analysis information:" << std::endl;
std::cout << " Name:  " << Name << std::endl;
std::cout << " Directory:  " << Directory << std::endl;

// Windows 
/*
Input_Dir            = Directory + "\\Input\\"   + Name + ".txt";
OutputMatlab_Dir     = Directory + "\\Output\\"  + Name + ".Matlab";
Info_Dir             = Directory + "\\Output\\"  + Name + ".inf";

FullFile_Dir         = Directory + "\\Output\\"  + Name + ".Res";
HistoryFile_Dir      = Directory + "\\Output\\"  + Name + ".His";
TransferFunction_Dir = Directory + "\\Output\\"  + Name + ".TRF";
*/

// Linux

Input_Dir            = Directory + "/input/"   + Name + ".txt";
OutputMatlab_Dir     = Directory + "/output/"  + Name + ".Matlab";
Info_Dir             = Directory + "/output/"  + Name + ".inf";

FullFile_Dir         = Directory + "/output/"  + Name + ".Res";
HistoryFile_Dir      = Directory + "/output/"  + Name + ".His";
TransferFunction_Dir = Directory + "/output/"  + Name + ".TRF";

// Final Directories
std::cout<< "  Input dir:   " << Input_Dir << std::endl;
std::cout<< "  Output dir:  " << Info_Dir  << std::endl;

Addressfile.close();

}




