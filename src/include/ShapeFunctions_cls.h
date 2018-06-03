


#ifndef SHAPEFUNCTIONS_CLS_H
#define SHAPEFUNCTIONS_CLS_H

namespace main_ns
{
  
namespace ShapeFunctions_ns
{

template<int NInt, int NNode>
class ShapeFunctions_cls(){

public: 

double XInt[NInt] ;    // ABSCISSAE of the integration points
double WInt[NInt] ;    // WEIGHTS of the integration points

double X1;            // Location of the integration 
double Fn[NNode];     // Shape function based on the number of nodes
double DFXI[NNode] ;      // The differential of the shape function

ShapeFunctions_cls();
void Retrieving_Gauss_Points_fn();

void ShapeFunc_1D_2N (double& x1, double& F0, double& F1);
void Dif_ShapeFunc_1D_2N ( double& DF0, double& DF1);
void ShapeFunc_1D_3N (double& x1, double& F0, double& F1, double& F2);
void Dif_ShapeFunc_1D_3N (double& x1, double& DF0, double& DF1, double& DF2);


  
}  
}

}

#endif
