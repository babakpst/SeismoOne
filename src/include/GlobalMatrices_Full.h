
//#include "ShapeFunctions.h"
//#include "Assemble.h"

#ifndef GLOBALMATRICES_FULL_H
#define GLOBALMATRICES_FULL_H

namespace main_ns
{

namespace solver_full_ns
{

class Matrices_full_cls{

public:
void MassDampStiffS_1D_first_Full ( int& iel, int& NEl, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT);
void MassDampStiffS_1D_second_Full ( int& iel, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT);
void GlobalMatrices_Full( int & NInt, int& NDim, int& NDOF, int& NNode, int& NEl, int *& MTel, int **& INod, int **& ID, double **& XYZ, double **& M, double **& C, double **& K );
void DRM_Matrices_Full( int & NNBndry, int &NNLayer, double **& K, double **& C, double **& M, double **& K_eb, double **& C_eb, double **& M_eb, int *&ND_e, int *&ND_b ) ;
}


}


}


#endif
