
#include "../include/Solver_Skyline.h"
using namespace std;



//***************************************************************************************************************************************************
// reduction of the effective stiffness matrix to a triangular one for Gauss-Elimination method using the skyline method
//***************************************************************************************************************************************************

void Skyline ( int & NEqM, int & NEl, int & NNode, int & NDOF, int *& NTK, int**& INod, int**& ID, int *& JD ) {

int i,j,k,l,iel;  // Loop indices
int NEqEl;

int    * ND;     // element constraints

// -----------------------------------------------------------------------------------------------------

NEqEl  = NDOF * NNode ;
ND = new int [NEqEl];

  for (i = 0; i < NEqEl; i++){
    ND[i]=0.0;
  }

  for (i=0; i < NEqM; i++){
    NTK[i] = i ;
  }

  for (iel=0; iel < NEl; iel++){ 
    for (i = 0; i<NNode; i++){
      k = INod[i][iel];
        for (j=0; j < NDOF; j++){
          ND[j * NNode + i]=ID[k][j] ;
        }
    }

    for (l =0; l<NEqEl; l++){
      for (k=0; k<NEqEl; k++){ 
        i=ND[l];
        j=ND[k];
        //if ((i==0) ||  (j==0)) continue;
        if (i>j) continue;
        if (i<NTK[j]) NTK[j]=i;
      }
    }
  }

JD[0]=0;

  for ( i=1; i<NEqM; i++){ 
    JD[i]=JD[i-1]+i+1-NTK[i];
  }

cout <<"End function skyline"<< endl; 

}

//***************************************************************************************************************************************************
// reduction of the effective stiffness matrix to a triangular one for Gauss-Elimination method using skyline method
//***************************************************************************************************************************************************

void Reduce_Skyline ( int& NEqM,      double *& K_S, int *& NTK, int *& JD, ofstream& info )
{

int N,N1,K,K1,I,I1,J,KJ,KK,KI,IJ;
double Fac;

N=NEqM;
N1=N-1;

  for (K=0;K<N1;K++) {
    K1=K+1;
      for (J=K1;J<=N;J++){
        if (K < (NTK[J])) continue;
        KJ=JD[J]+K-J;
        KK=JD[K];
        I1=K1;
        if (K1 < (NTK[J])) I1=NTK[J];
          for (I=I1;I<=J;I++){
            if (K > (NTK[I])) {
              KI=JD[I]+K-I;
              Fac=K_S[KI]/K_S[KK];
              IJ=JD[J]+I-J;
              info << "i " << I << " j " << J<< " Fac " << Fac << " K " << K_S[IJ] << " IJ " << IJ << endl;
              K_S[IJ]=K_S[IJ]-K_S[KJ]*Fac;
            }
          }
      }
  }


/*

n=NEqM;
n1=n-1;


  for (k=1;k<=n1;k++) {
    k1 = k + 1 ;
      for (j=k1;j<=n;j++){
        if (k < (NTK[j-1]+1)) continue;
        kj=JD[j-1]+k-j+1;
        kk=JD[k-1]+1;
        i1=k1;
        if (k1 < (NTK[j-1]+1)) i1=NTK[j-1]+1;
          for (i=i1;i<j;i++){
cout<< k << " " << (NTK[i-1]+1) << endl;
            if (k > (NTK[i-1]+1)) {
              ki=JD[i-1]+k-i+1;
              Fac=K_S[ki-1]/K_S[kk-1];
              ij=JD[j-1]+i-j+1;
              Check << "i " << i << " j " << j<< " Fac " << Fac << " K " << K_S[ij-1] << " IJ " << ij << endl;
              cout << "i " << i << " j " << j<< " Fac " << Fac << " K " << K_S[ij-1] << " IJ " << ij << endl;
              K_S[ij-1]=K_S[ij-1]-K_S[kj-1]*Fac;
            }
          }
      }
  }



*/

}


//***************************************************************************************************************************************************
// Solver: Gauss-Elimination
//***************************************************************************************************************************************************
void Gauss_El_Skyline ( int *& NTK, int *& JD, int& NEqM, double *& UN, double *& K_S)
{


int k, j, KJ, i, N, N1, KK, K1, NN;          // temporary variables

N = NEqM ;
N1 = N - 1 ;

  for (k=0;k<N1; k++) {
    K1 = k + 1 ;
      for (int j=K1;j<N;j++) {
        if (k < NTK[j]) continue;
          KJ=JD[j]+k-j;
          KK=JD[k];
          UN[j]=UN[j]-UN[k]*K_S[KJ]/K_S[KK];
      }
  }

NN=JD[N];

UN[N]=UN[N]/K_S[NN];

  for (i=0;i<N1;i++) {
    k = N - i ;
    K1 = k + 1 ;
      for (j=K1;j<N;j++) {
        if ( k >= NTK[j]) {
          KJ=JD[j]+k-j;
          UN[k]=UN[k]-K_S[KJ]*UN[j];
        }
      }
    KK=JD[k];
    UN[k]=UN[k]/K_S[KK];
  }

}




void Matrix_Multiplication ( int *& NTK, int *& JD, double * &M1, double * &M2, double * &M3, int NEqM ){

int I,J,K,II,JJ,IJ;

  for (II=0; II<NEqM; II++){
      for (JJ=0;JJ<NEqM;JJ++){
        I=II;
        J=JJ;
          if (I>J) {
            K=I;
            I=J;
            J=K;
          }
        if (I<NTK[J]) continue;
        IJ=JD[J]+I-J;
        M3[II] += M1[IJ]*M2[JJ];
      }
  }

}

//***************************************************************************************************************************************************
// Newmark's algorithm for marching in time
//***************************************************************************************************************************************************
void Newmark_Skyline ( double & L, int & Wave_Type, int & Wave_Func, int &NStep, int& NEqM, int& LoadType, double& Gama, double& Beta, double& DT, double& Alpha, double *& M_S, double *& C_S, double *& K_S, double *& F, double **& PMat, double **& XYZ, ofstream& FullSol, ofstream& History, int *&ND_e, int *&ND_b, int *&Nodal_History, int *&JD, int *& NTK, ofstream& Check ) { 

int    i, ij ;                  // Loop indices

double A0, A1, A2, A3, A4, A5;  // Newmark variables
double Time;                    // time
double Total_Time;              // time
double LoadFactor;              // Load Factor at each time step
double E;                       // Elastic Modulus
double Rho;                     // density
double c ; 
double Initial_Time ; // Starting time of the simulation

double * UN;                    // temporay arrays
double * U;                     // temporay arrays
double * UD;                    // temporay arrays
double * UDD;                   // temporay arrays
double * Temp;                  // temporay arrays

bool InitialTime = false;


cout << "Newmark solver -- skyline" << endl;

E   = PMat[0][0];
Rho = PMat[0][1];
c = sqrt (E/Rho);

UN   = new double [NEqM];
U    = new double [NEqM];
UD   = new double [NEqM];
UDD  = new double [NEqM];
Temp = new double [NEqM];

// Newmark constants
A0 = 1.0 / ( Beta * DT * DT ) ;
A1 = Gama   / ( Beta * DT ) ;
A2 = 1.0 / ( Beta * DT ) ;
A3 = 1.0 / ( 2.0 * Beta ) - 1.0 ;
A4 = Gama   / Beta - 1.0 ;
A5 = DT * ( Gama / ( 2.0 * Beta ) - 1.0 ) ;

// Initializing displacement, velocity and acceleration
  for (i=0;i<NEqM;i++) {
    UN[i]  = 0.0 ;
    U[i]   = 0.0 ;
    UD[i]  = 0.0 ;
    UDD[i] = 0.0 ;
  }

// Effective stiffness matrix 
std::cout << "Effective matrix ..." << endl;

  for ( i=0; i<JD[NEqM-1]; i++) {
    K_S[i] = K_S[i] + A0 * M_S[i] + A1 * C_S[i];
  }


std::cout << "Reduce effective matrix ..." << endl;
// Reduction the coefficient matrix ( Effective Stiffness Matrix )
Reduce_Skyline ( NEqM,      K_S, NTK, JD, Check );



//  for (k=0; k<NEqM; k++){  
//    cout<<JD[k]<<"   "<< NTK[k]<< endl;
//  }
/*
cout<< "K22 " << K_S[2] <<"K3 " << K_S[3] <<"K4 " << K_S[4] <<"K5 " << K_S[5] << endl;

  for (k=0; k<10; k++){
    cout<< "This is I :" << k << endl;
    i=k;
      for ( l=0; l<10; l++){
        j=l;
          if (i>j) {
            ij = i;
            i  = j;        
            j  = ij;
          }
          if ( j<NTK[j]) TE =0.0;
          else {
            ij = JD[j]+i-j;
            TE = K_S[ij] ;
          }
        Check << TE<< setw(20);
      }
    Check << endl ;
  }

cin.get();
*/

Initial_Time = -L / c; 
std::cout << "Initial time = " << Initial_Time << endl;

  // Solve the PDE for each step -----------------------------------------------
  for (int IStep=0;IStep<NStep+1;IStep++) {

    Total_Time = IStep * DT; 
    Time  = Initial_Time + IStep * DT ;  // Time STEP
    std::cout << "Time Step:  " << IStep << "  Time: " << Time << "  Total time:" << Total_Time << endl;

    // Update displacements, velocity and acceleration
      for (i=0; i<NEqM; i++) {
        Temp[i] = UD[i] ;
        U[i]    = UN[i] - U[i] ;
        UD[i]   = A1 * U[i] - A4 * Temp[i] - A5 * UDD[i] ;
        UDD[i]  = A0 * U[i] - A2 * Temp[i] - A3 * UDD[i] ;
        U[i]    = UN[i] ;
        UN[i]   = 0.0 ; 
      }

    /*
    // Recording the full results
    FullSol << "Time:  " << Time << "\n";
    HistorySolution ( NJ, Time, Alpha, P, E, Rho, A,   U_EX, XYZ) ;
      // full results
      for (ij=0;ij<NJ;ij++) {
        if ( ID [ij][0] < 0 ) u = 0.0;
        else u = U[ij];
        FullSol << setw(15) << XYZ[ij][0] << setw(15) << u << setw(15) << U_EX[ij] << "\n";
      }
    */

      // Effective force - stored in UN
      for (i=0;i<NEqM; i++) { // find the coefficient of the M matrix
        Temp[i] = A0 * U[i] + A2 * UD[i] + A3 * UDD[i] ;
      }



    Matrix_Multiplication ( NTK, JD, M_S, Temp, UN, NEqM );
/*
      for (i=0;i<NEqM;i++) { // multiplying the vector by the mass matrix
        TE = 0.0;
          for (int j=0;j<NEqM;j++) {
            TE += M[i][j] * Temp[j];
          }
        UN[i] = TE;
      }
*/



      for (i=0;i<NEqM; i++) {
        Temp[i] = A1 * U[i] + A4 * UD[i] + A5 * UDD[i] ;
      }

    Matrix_Multiplication ( NTK, JD, C_S, Temp, UN, NEqM );

/*
      for (i=0;i<NEqM;i++) {
        TE = 0.0;
          for (int j=0;j<NEqM;j++) {
            TE  += C[i][j] * Temp[j] ;
          }
        UN[i] +=  TE;
      }
*/


    // Adding load at this time step
      if ( LoadType == 0 ) {
        LoadFactor = LoadFunction ( Time, Alpha, P);    // Pressure load
        for (ij=0; ij<NJ; ij++) {
          UN[ij] = UN[ij] - F[ij] * LoadFactor;
        }
      }
      else if ( LoadType == 1 ) { //        DRM_Load ();
        F[0] = 0;
 
        DRM_Loads_Implicit ( alpha1, alpha2, Time, NDim, NNBndry, NNLayer,              Wave_Type, Wave_Func, amplitude, c, UN, XYZ, NoBndry_DRM, NoLayer_DRM, M_eb, C_eb, K_eb, ND_e, ND_b ); 
      }

    // Check whether the initial time is small enough
      if ( IStep == 0) {
        for (i=0; i<NEqM; i++ ) {
          if (UN[i]!= 0) {
            InitialTime = true;
            break;
          }
        }
        if (InitialTime == true) {
          cout << "WARNING: REDUCE THE INITIAL TIME" << endl;
          return;
        }
      }

    // SOLVE
    Gauss_El_Skyline ( NTK, JD, NEqM, UN, K_S);

    // time history of the solution at some particular nodes
    History << setw(6) << Total_Time;
      for (i=0;i<Dis_History;i++) {
        History << setw(20) << UN[Nodal_History[i]]  ;
      }
    History  << endl;


    //FullSol << "Time= " << Time << endl;
      for (i = 0; i<NEqM; i++){
        FullSol << UN[i]<< setw(20);
      }
    FullSol << endl;

  }

}




