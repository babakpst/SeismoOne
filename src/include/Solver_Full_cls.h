
#include <iostream>
#include <fstream>

#include "../include/solver_cls.h"
#include "../include/Discretization_cls.h"

#ifndef SOLVER_FULL_H
#define SOLVER_FULL_H

namespace main_ns{

namespace solver_full_ns{
  
class solver_full_cls: public main_ns::solver_ns::solver_cls{


public:

	explicit solver_full_cls(main_ns::discretization_ns::discretization_cls*,main_ns::model_ns::model_cls*);


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

