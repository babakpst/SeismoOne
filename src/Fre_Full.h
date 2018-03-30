
#include "Variables.h"
#include "Load.h"

#ifndef FRE_FULL
#define FRE_FULL

void LDLT_Freq ( int& NEqM, double **& K_Eff);
void Substitute_Freq ( int& NEqM, double *& RHS, double **& K_Eff);
void Transfer_Full ( double & alpha1, double & alpha2, double & L, int & Wave_Type, int& NEqM, double **& M, double **& C, double **& K, double *& F, double **& PMat, double **& XYZ, int *&ND_e, int *&ND_b, int *&Nodal_History, ofstream& info, ofstream& TransferFunc );

#endif

