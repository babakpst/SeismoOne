


#ifndef SHAPEFUNCTIONS_CLS_H
#define SHAPEFUNCTIONS_CLS_H

namespace main_ns
{
  
namespace ShapeFunctions_ns
{

class ShapeFunctions_cls(){

public: 
ShapeFunctions_cls();

void ShapeFunc_1D_2N (double& x1, double& F0, double& F1);
void Dif_ShapeFunc_1D_2N ( double& DF0, double& DF1);
void ShapeFunc_1D_3N (double& x1, double& F0, double& F1, double& F2);
void Dif_ShapeFunc_1D_3N (double& x1, double& DF0, double& DF1, double& DF2);
void GAUSS_Quad_POINTS ( int& NInt, Gauss& GAUSS_POINTS);
  
}  
}

}

#endif
