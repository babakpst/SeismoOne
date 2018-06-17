
#include <iostream>
#include <fstream>

#include "../include/create_global_matrices_cls.h"
#include "../include/Discretization_cls.h"
#include "../include/ShapeFunctions_cls.h"
#include "../include/ShapeFunctions_FirstOrder_cls.h"
#include "../include/ShapeFunctions_SecondOrder_cls.h"

#ifndef CREATE_FULL_MATRICES_CLS_H
#define CREATE_FULL_MATRICES_CLS_H

namespace main_ns{

namespace Matrices_Full_ns{
  
class Matrices_Full_cls: public main_ns::Matrices_ns::Matrices_cls{

  int MType;              // Material type
  double E;               // elastic modulus
  double Rho;             // density

protected:
	void allocating_global_matrices_fn();
  void allocating_local_matrices_fn();
  

public:
	Matrices_Full_cls( main_ns::discretization_ns::discretization_cls*, 
	                            main_ns::model_ns::model_cls*);

  virtual void assembling_local_matrices_into_global_matrices_fn();  

  void compute_elemental_matrices_fn(const int*, const double*, const double*);


/*
//void Reduce_Full (int& NEqM, double **& K, ofstream& Check);
void LDLT ( int& NEqM, double **& K); 
void Substitute ( int& NEqM, double *& UN, double **& K);
//void Gaussian ( int& NEqM, double *& UN, double **& K);
void Newmark_Full ( double & L, int & Wave_Type, int & Wave_Func, int &NStep, int& NEqM, 
                    int& LoadType, double& Gama, double& Beta, double& DT, double& Alpha, 
                    double **& M, double **& C, double **& K, double *& F, double **& PMat, 
                    double **& XYZ, ofstream& FullSol, ofstream& History, int *&ND_e, int *&ND_b, 
                    int *&Nodal_History );
*/

};
}
}

#endif

