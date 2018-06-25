


#ifndef LOAD_H
#define LOAD_H

namespace main_ns
{

namespace Solver_ns
{

class apply_seismic_loads_to_the_domain
{

double LoadFunction ( double& Time, double& Alpha, double& P) ;
void HistorySolution ( int& NJ, double& TIME, double& Alpha, double& P, double& E, double& Rho, double& A,   double *& U_EX, double **& XYZ);
void DRM_PointValues (int & Wave_Func, double & amplitude, double & Time, double & x, double & c, double & omega, double & alpha1, double & alpha2, double & u, double & v, double & a);
void DRM_Loads_Implicit ( double & alpha1, double & alpha2, double & Time, int NDim, int NNBndry, int NNLayer, int& Wave_Type, int& Wave_Func, double & amplitude, double & c, double *& UN, double **& XYZ, int *& NoBndry_DRM, int *& NoLayer_DRM, double **& M_eb, double **& C_eb, double **& K_eb, int *&ND_e, int *&ND_b ); 
void DRM_PointValues_Freq ( double & amplitude, double & x, double & c, double & omega, double & u_R, double & u_I );

} // load class
} // Solver_ns
} // main_ns

#endif
