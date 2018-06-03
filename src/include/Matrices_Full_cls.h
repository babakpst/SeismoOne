
#include <iostream>
#include <fstream>

#include "../include/Matrices_cls.h"
#include "../include/Discretization_cls.h"
#include "../include/ShapeFunctions_cls.h"

#ifndef SOLVER_FULL_H
#define SOLVER_FULL_H

namespace main_ns{

namespace Matrices_Full_ns{
  
class Matrices_Full_cls: public main_ns::Matrices_ns::Matrices_cls{

  void allocating_local_matrices_fn();
	void allocating_global_matrices_fn();
  void compute_elemental_matrices_fn();

public:
	explicit Matrices_Full_cls(main_ns::discretization_ns::discretization_cls*, 
	                           main_ns::model_ns::model_cls*);
  void assembling_local_matrices_into_global_matrices_fn();



/*
//void Reduce_Full (int& NEqM, double **& K, ofstream& Check);
void LDLT ( int& NEqM, double **& K); 
void Substitute ( int& NEqM, double *& UN, double **& K);
//void Gaussian ( int& NEqM, double *& UN, double **& K);
void Newmark_Full ( double & L, int & Wave_Type, int & Wave_Func, int &NStep, int& NEqM, int& LoadType, double& Gama, double& Beta, double& DT, double& Alpha, double **& M, double **& C, double **& K, double *& F, double **& PMat, double **& XYZ, ofstream& FullSol, ofstream& History, int *&ND_e, int *&ND_b, int *&Nodal_History );
*/

};
}
}

#endif

