
// include files

using namespace std;

#include "../include/Address.h"


// Constructor
main_ns::address_ns::address_cls::address_cls(){
std::cout << " --------------- SeismoOne ---------------" << endl;
std::cout << " ----- Developed by: Babak Poursartip ----" << endl;
std::cout << " -------------- version 3.0 --------------" << endl;
std::cout << endl;
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
std::cout << " opening the address file ..." << endl;
ifstream Addressfile;
Addressfile.open ("Address.txt", ios::in );

// Reading simulation information
std::cout << " Reading file information ..." << endl;
getline (Addressfile,TempS);
Addressfile >> Name ; // Input file name

getline (Addressfile,TempS);
getline (Addressfile,TempS);
Addressfile >> Directory ; // Directory of the input file

std::cout << endl;
std::cout << " Analysis information:" << endl;
std::cout << " Name:  " << Name << endl;
std::cout << " Directory:  " << Directory << endl;

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
std::cout<< "  Input dir:   " << Input_Dir << endl;
std::cout<< "  Output dir:  " << Info_Dir  << endl;

Addressfile.close();

}




