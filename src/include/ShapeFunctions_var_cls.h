


#ifndef SHAPEFUNCTIONS_CLS_H
#define SHAPEFUNCTIONS_CLS_H

namespace main_ns
{
  
namespace ShapeFunctions_ns
{

template<int NInt, int NNode>
class ShapeFunctions_variables_cls(){

public: 

double X1;            // Location of the integration 
double Fn[NNode];     // Shape function based on the number of nodes
double DFXI[NNode];  // The differential of the shape function

double XInt[NInt] ;    // ABSCISSAE of the integration points
double WInt[NInt] ;    // WEIGHTS of the integration points

ShapeFunctions_variables_cls();
~ShapeFunctions_variables_cls();
  
}  
}

}

#endif
