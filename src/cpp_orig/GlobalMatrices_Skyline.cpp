

//***************************************************************************************************************************************************
// Copy the submatrices for DRM loads.  
//***************************************************************************************************************************************************
void DRM_Matrices_Skyline( int & NNBndry, int &NNLayer, double *& K_S, double *& C_S, double *& M_S, double **& K_eb, double **& C_eb, double **& M_eb, int *&ND_e, int *&ND_b , int *& JD)  
{

int ij,i,j,l,n;   // Loop indices

// - Code ---------------------------------------------------------------------
std::cout << "Create DRM matrices ..." << endl; 

  for ( l = 0; l < NNLayer * NDim; l++) {
    for ( n = 0; n < NNBndry * NDim; n++) {

      i = ND_e [ l ] ;
      j = ND_b [ n ] ;


      if (i>j) {
        ij = i;
        i  = j;
        j  = ij; 
        std::cout << " i and j replaced"<< endl;
      }

      ij = JD[j]+i-j;

//cout<< "DRM" << ij << " K_S  "<< K_S[ ij ] <<  " M_S "<< M_S[ ij ] << " C_S  "<< C_S[ ij ] <<  endl;
   
      K_eb[l][n] = K_S[ ij ];
      C_eb[l][n] = C_S[ ij ];
      M_eb[l][n] = M_S[ ij ];

//      cin.get();
    }
  }

std::cout << "DRM matrices created." << endl; 


}








