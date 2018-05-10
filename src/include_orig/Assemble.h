#include "Variables.h"

#ifndef ASSEMBLE
#define ASSEMBLE

void AssembleMassDampingStiffForceFull (int& NEqEl, int *& ND, double **& Ke, double **& Ce, double **& Me, double **& K, double **& C, double **& M);
void AssembleMassDampingStiffForceSkyline (int& NEqEl, int *& ND, double **& Ke, double **& Ce, double **& Me, double *& K_S, double *& C_S, double *& M_S, int *& JD);

#endif
