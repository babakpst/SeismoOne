
#include "../include/Assemble.h"


//***************************************************************************************************************************************************
// Assembling element matrices in the global matrices
//***************************************************************************************************************************************************
void AssembleMassDampingStiffForceFull  (int& NEqEl, int *& ND, double **& Ke, double **& Ce, double **& Me, double **& K, double **& C, double **& M) 
{

int L;
int N;

  for (int i=0;i<NEqEl;i++) {
    for (int j=0;j<NEqEl;j++) {
      L = ND [ i ] ;
      N = ND [ j ] ;
        if ( L == -1 || N == -1 ) continue;
      K [L][N] = K [L][N] + Ke [i][j];
      C [L][N] = C [L][N] + Ce [i][j];
      M [L][N] = M [L][N] + Me [i][j];
    }
  }

// std::cout << "C  " <<  C [L][N] << "Ce" << Ce [i][j] << endl;
//cin.get();

}



//***************************************************************************************************************************************************
// Assembling element matrices in the global matrices
//***************************************************************************************************************************************************
void AssembleMassDampingStiffForceSkyline  (int& NEqEl, int *& ND, double **& Ke, double **& Ce, double **& Me, double *& K_S, double *& C_S, double *& M_S, int *& JD ) 
{

int l,n, i, j,ij;

  for ( l=0;l<NEqEl;l++) {
    for ( n=0;n<NEqEl;n++) {
      i = ND [ l ] ;
      j = ND [ n ] ;
        if (i>j) continue;
        ij=JD[j]+i-j; 

        K_S[ij]=K_S[ij]+Ke[l][n]; 
        C_S[ij]=C_S[ij]+Ce[l][n]; 
        M_S[ij]=M_S[ij]+Me[l][n]; 
    }
  }

}
