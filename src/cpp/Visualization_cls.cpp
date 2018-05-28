
#include "../include/Visualization.h"


main_ns::visualization_ns::visualization_cls::visualization_cls(
                               const main_ns::address_ns::address_cls* a,
                               const main_ns::model_ns::model_cls* b,
                               const main_ns::discretization_ns::discretization_cls* c):
                               address (a),
                               model(b),
                               discretized_model(c)
                               {
std::cout << " Writing the input file for Matlab visualizer script ..." << std::endl;

//address=a;
//model=b;
//discretized_model=c;

// Output file for Matlab for visualization
std::cout << " Opening the results file for Matlab ..."<< std::endl;
OutputMatlab.open (address->OutputMatlab_Dir.c_str(), std::ios::out);

}

/*
###################################################################################################
Purpose: This function writes down all the required information for the Matlab script to 
         visualize the results.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin
================================= V E R S I O N ===================================================
V0.00: 05/20/2018 - Subroutine initiated.
V0.01: 05/21/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/


void main_ns::visualization_ns::visualization_cls::MatlabOutput_fn(){

OutputMatlab << model->NMat << std::endl;
OutputMatlab << model->NEl << std::endl;
OutputMatlab << model->Dis_History << std::endl;
OutputMatlab << model->DT << std::endl;
OutputMatlab << model->NStep << std::endl;
OutputMatlab << model->L << std::endl;
OutputMatlab << discretized_model->NEqM << std::endl;

  for (int i = 0; i<model->NMat; i++){
      for (int j = 0; j<model->NPM; j++){
        OutputMatlab << std::setw(20) << model->PMat[i][j] ;
      }
    OutputMatlab << std::endl;
  }

  for (int i = 0; i<model->NMat; i++){
    OutputMatlab << std::setw(20) << model->Length[i] ;
  }
OutputMatlab << std::endl;

  for (int i = 0; i<model->Dis_History; i++){
    OutputMatlab << std::setw(20) << model->Loc_History [i]  ;
  }
OutputMatlab << std::endl;

  for (int i = 0; i<model->NJ; i++){
    OutputMatlab << discretized_model->XYZ[i][0]  << std::endl;
  }

OutputMatlab << std::endl;

std::cout << " Matlab input file created successfully." << std::endl;
OutputMatlab.close();

}
