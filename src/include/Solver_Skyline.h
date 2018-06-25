
#include "Variables.h"
#include "Load.h"

#ifndef SOLVER_SKYLINE
#define SOLVER_SKYLINE

void Skyline ( int & NEqM, int & NEl, int & NNode, int & NDOF, int *& NTK, int**& INod, int**& ID, int *& JD ) ;
void Reduce_Skyline ( int& NEqM,      double *& K_S, int *& NTK, int *& JD, ofstream& info );
void Gauss_El_Skyline ( int *& NTK, int *& JD, int& NEqM, double *& UN, double *& K_S) ;
void Matrix_Multiplication ( int *& NTK, int *& JD, double * &M1, double * &M2, double * &M3, int NEqM );


void Newmark_Skyline ( double & L, int & Wave_Type, int & Wave_Func, int &NStep, int& NEqM, int& LoadType, double& Gama, double& Beta, double& DT, double& Alpha, double *& M_S, double *& C_S, double *& K_S, double *& F, double **& PMat, double **& XYZ, ofstream& FullSol, ofstream& History, int *&ND_e, int *&ND_b, int *&Nodal_History, int *&JD, int *& NTK, ofstream& Check );

#endif



