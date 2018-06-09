

#include <math.h> 
#include <iostream>


#ifndef SHAPEFUNCTIONS_CLS_H
#define SHAPEFUNCTIONS_CLS_H

namespace main_ns
{
  
namespace ShapeFunctions_ns
{

class ShapeFunctions_cls{

public: 

int NInt;         // Number of Integration points (quadrature)
int NNode;        // Number of nodes in the element

double x1;       // Location of the integration 
double* Fn;      // Shape function based on the number of nodes
double* DFXI;    // The differential of the shape function

double* XInt;    // ABSCISSAE of the integration points
double* WInt;    // WEIGHTS of the integration points

ShapeFunctions_cls(int, int);

void Retrieving_Gauss_Points_fn();

virtual void ShapeFunctions (void) = 0;
virtual void DifferentialOfShapeFunctions (void)  = 0;

virtual ~ShapeFunctions_cls();
  
}; 
}
}

#endif




  
