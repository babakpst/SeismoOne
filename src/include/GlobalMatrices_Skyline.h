
/*
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
*/

#include "Variables.h"
#include "ShapeFunctions.h"
#include "Assemble.h"


#ifndef GLOBALMATRICES_SKYLINE_H
#define GLOBALMATRICES_SKYLINE_H

void MassDampStiffS_1D_first_Skyline ( int& iel, int& NEl, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT);
void MassDampStiffS_1D_second_Skyline ( int& iel, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT);
void GlobalMatrices_Skyline ( int & NInt, int& NDim, int& NDOF, int& NNode, int *& MTel, int **& INod, int **& ID, double **& XYZ, double *& M_S, double *& C_S, double *& K_S, int *& JD) ;
void DRM_Matrices_Skyline( int & NNBndry, int &NNLayer, double *& K_S, double *& C_S, double *& M_S, double **& K_eb, double **& C_eb, double **& M_eb, int *&ND_e, int *&ND_b , int *& JD)  ;

#endif
