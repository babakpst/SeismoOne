
#include <iostream>
#include <fstream>

#include "../include/create_global_matrices_cls.h"
#include "../include/Discretization_cls.h"

#include "../include/ShapeFunctions_cls.h"
#include "../include/ShapeFunctions_FirstOrder_cls.h"
#include "../include/ShapeFunctions_SecondOrder_cls.h"

//#include "../include/assemble_local_to_global.h"

#ifndef CREATE_SKYLINE_MATRICES_CLS_H
#define CREATE_SKYLINE_MATRICES_CLS_H

namespace main_ns{

namespace Matrices_ns{
  
class Matrices_Skyline_cls: public main_ns::Matrices_ns::Matrices_cls, 
                         public main_ns::Matrices_ns::assemble_local_to_global_cls{

protected:

    int * JD;             // Skyline matrix
    int * NTK;            // Skyline matrix

    double * K_S;         // global stiffness matrix -skyline
    double * C_S;         // global damping matrix -skyline
    double * M_S;         // global mass matrix -skyline

	virtual void allocating_global_matrices_fn();
  

public:
	Matrices_Full_cls(main_ns::discretization_ns::discretization_cls*,main_ns::model_ns::model_cls*);
  
  virtual void assemble_local_to_global_fn(void)=0;

};
}
}

#endif

