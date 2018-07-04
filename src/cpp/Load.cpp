

#include "../include/Load.h"

main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls::
    apply_seismic_loads_to_the_domain_cls(){};

/*
###################################################################################################
Purpose: This function computes the load factor to apply the pressure load on the surface.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/26/2018 - Subroutine initiated.
V1.00: 06/26/2018 - Compiled successfully.

###################################################################################################
*/

double main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls::
    LoadFunction(const double Time, const double Alpha, const double P)
{
  if (Time < 2.0 * pi / Alpha)
    LoadFactor = -P * sin(Alpha * Time);
  else
    LoadFactor = 0.0;

  return (LoadFactor);
}










///////////////////////////////////////////////////////
/*
###################################################################################################
Purpose: This function compuates the DRM point loads at each load.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/26/2018 - Subroutine initiated.
V1.00: 06/27/2018 - Compiled successfully.

###################################################################################################
*/

void main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls::
    DRM_PointValues(PointLoad Load)
{

  int TotalCycle;
  int Direction;

  double arg1; // wave phase for the moving in the positive direction
  double arg2; // wave phase for the moving in the negative direction

  double LowerLimit; // Lower limit of the phase
  double UpperLimit; // Upper limit of the phase

  double wr;    // Central frequency of the Ricker pulse
  double t_max; // Loading duration
  double fr;

  double ur;

  // The analytical solution is u (x,t) = Ui (f_inc (t - x/c) + f_inc (t + x/c) )

  // phases

  switch (Load.Wave_Func)
  {
  case 0: // sine function

    wr = 2.0 * pi * Load.omega; // characteristic central circular frequency

    arg1 = Load.Time - Load.x / Load.c; // positive direction phase - incident wave
    arg2 = Load.Time + Load.x / Load.c; // negative direction phase - reflected wave

    arg1 *= wr;
    arg2 *= wr;

    LowerLimit = Load.alpha1;
    UpperLimit = Load.alpha2;

    if (LowerLimit <= arg1 && arg1 <= UpperLimit)
    {
      Load.u +=           Load.amplitude * sin(arg1);
      Load.v += wr *      Load.amplitude * cos(arg1);
      Load.a -= wr * wr * Load.amplitude * sin(arg1);
    }

    if (LowerLimit <= arg2 && arg2 <= UpperLimit)
    {
      Load.u +=           Load.amplitude * sin(arg2);
      Load.v += wr *      Load.amplitude * cos(arg2);
      Load.a -= wr * wr * Load.amplitude * sin(arg2);
    }

    break;
  case 1: // Ricker

    TotalCycle = (int)(Load.alpha1);
    Direction  = (int)(Load.alpha2);

    fr = Load.omega;    // Central frequency of Ricker pulse
    wr = 2.0 * pi * fr; // characteristic central circular frequency

    //arg1 = wr * Time - wr * x/c ; // positive direction phase - incident wave
    //arg2 = wr * Time + wr * x/c ; // negative direction phase - reflected wave

    arg1 = Load.Time - Load.x / Load.c; // positive direction phase - incident wave
    arg2 = Load.Time + Load.x / Load.c; // negative direction phase - reflected wave

    //wr = omega ;           // characteristic central circular frequency

    t_max = 6.0 * sqrt(6.0) / wr; // duration of loading

    LowerLimit = 0.0;
    UpperLimit = t_max;

    if (LowerLimit <= arg1 && arg1 <= UpperLimit)
    {

      ur = wr * arg1 - 3.0 * sqrt(6.0);

      Load.u += Load.amplitude * ((0.25 * ur * ur - 0.5) * exp(-0.25 * ur * ur) - 13.0 * exp(-13.5)) / (0.5 + 13.0 * exp(-13.5));
      Load.v += Load.amplitude * (wr * (0.75 * ur - 0.125 * pow(ur, 3.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
      Load.a += Load.amplitude * (pow(wr, 2.0) * (0.75 - 0.75 * pow(ur, 2.0) + 0.0625 * pow(ur, 4.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
    }
    else
    {
      for (int j = 2; j <= TotalCycle; j++)
      {
        arg1 = Load.Time - Load.x / Load.c - (j - 1) * t_max;
        if (LowerLimit <= arg1 && arg1 <= UpperLimit)
        {
          ur = wr * arg1 - 3.0 * sqrt(6.0);

          Load.u += (pow(Direction, (j - 1))) * Load.amplitude * ((0.25 * ur * ur - 0.5) * exp(-0.25 * ur * ur) - 13.0 * exp(-13.5)) / (0.5 + 13.0 * exp(-13.5));
          Load.v += (pow(Direction, (j - 1))) * Load.amplitude * (wr * (0.75 * ur - 0.125 * pow(ur, 3.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
          Load.a += (pow(Direction, (j - 1))) * Load.amplitude * (pow(wr, 2.0) * (0.75 - 0.75 * pow(ur, 2.0) + 0.0625 * pow(ur, 4.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
        }
      }
    }

    if (LowerLimit <= arg2 && arg2 <= UpperLimit)
    {

      ur = wr * arg2 - 3.0 * sqrt(6.0);

      Load.u += Load.amplitude * ((0.25 * ur * ur - 0.5) * exp(-0.25 * ur * ur) - 13.0 * exp(-13.5)) / (0.5 + 13.0 * exp(-13.5));
      Load.v += Load.amplitude * (wr * (0.75 * ur - 0.125 * pow(ur, 3.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
      Load.a += Load.amplitude * (pow(wr, 2.0) * (0.75 - 0.75 * pow(ur, 2.0) + 0.0625 * pow(ur, 4.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
    }
    else
    {
      for (int j = 2; j <= TotalCycle; j++)
      {
        arg2 = Load.Time + Load.x / Load.c - (j - 1) * t_max;
        if (LowerLimit <= arg2 && arg2 <= UpperLimit)
        {
          ur = wr * arg2 - 3.0 * sqrt(6.0);

          Load.u += (pow(Direction, (j - 1))) * Load.amplitude * ((0.25 * ur * ur - 0.5) * exp(-0.25 * ur * ur) - 13.0 * exp(-13.5)) / (0.5 + 13.0 * exp(-13.5));
          Load.v += (pow(Direction, (j - 1))) * Load.amplitude * (wr * (0.75 * ur - 0.125 * pow(ur, 3.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
          Load.a += (pow(Direction, (j - 1))) * Load.amplitude * (pow(wr, 2.0) * (0.75 - 0.75 * pow(ur, 2.0) + 0.0625 * pow(ur, 4.0)) * exp(-0.25 * ur * ur)) / (0.5 + 13.0 * exp(-13.5));
        }
      }
    }

    break;
  default:
    std::cout << "Wave function is not defined." << std::endl;
  }
}

/*
###################################################################################################
Purpose: This function computes the load factor to apply the pressure load on the surface.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 07/04/2018 - Subroutine initiated.
V1.00: 07/04/2018 - Compiled successfully.

###################################################################################################
*/
void DRM_Loads_Implicit(const main_ns::Solver_ns::InputLoad* LoadPackage, 
                        const main_ns::Solver_ns::PointLoad* Load)
{

  // = Local Variables ================================================================================================================================
  double x; // The coordinate
  double u; // Analytical displacement
  double v; // Analytical velocity
  double a; // Analytical acceleration

  double *U_b;   // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.
  double *Ud_b;  // The vector that holds the analytical velocity at the boundary nodes. In this 1D problem, there is only one node. (Not needed because we do not have damping in the system)
  double *Udd_b; // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.
  double *F_b;   // The vector that holds the loads for the boundary nodes.

  double *U_e;   // The vector that holds the analytical displacement at the boundary nodes. In this 1D problem, there is only one node.
  double *Ud_e;  // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.  (Not needed because we do not have damping in the system)
  double *Udd_e; // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.
  double *F_e;   // The vector that holds the loads for the layer nodes.

  // = Function =======================================================================================================================================
  
  // Defining the required vectors
  U_b   = new double[LoadPackage.NNBndry * LoadPackage.NDim];
  Ud_b  = new double[LoadPackage.NNBndry * LoadPackage.NDim];
  Udd_b = new double[LoadPackage.NNBndry * LoadPackage.NDim];

  F_e   = new double[LoadPackage.NNLayer * LoadPackage.NDim];

  for (i = 0; i < LoadPackage.NNLayer * LoadPackage.NDim; i++)
  {
    F_e[i] = 0.0;
  }

  // Loop on the nodes on the DRM boundary to find out the analytical solution (In this case only one node)
  for (int i = 0; i < NNBndry; i++)
  {
    for (int j = 0; j < NDim; j++)
    {
      x = LoadPackage.XYZ[ LoadPackage.NoBndry_DRM[i]][j]; // Coordinate of the node
      u= v= a= 0.0;            // Initialize the values
      
      // Computing the analytical solution at this particular node
      // Remark: the one-dimensional wave-motion is identical for both SV and P waves.

      // The point loads are identical for both shear and pressure waves. 
      // Just for the sake of clarity, we wrtie it as follows.
      if (LoadPackage.Wave_Type == 0) // SV wave
        DRM_PointValues(Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); 
      else if (LoadPackage.Wave_Type == 1) // P wave
        DRM_PointValues(Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); 


      // Filling the analytical solution vector
      U_b  [i * LoadPackage.NNBndry * LoadPackage.NDim + j] = u;
      Ud_b [i * LoadPackage.NNBndry * LoadPackage.NDim + j] = v;
      Udd_b[i * LoadPackage.NNBndry * LoadPackage.NDim + j] = a;
    }
  }

  // Multiply the Mass, Damp, and Stiffness matrix by the vector
  for (int i = 0; i < LoadPackage.NNLayer * LoadPackage.NDim; i++)
  {
    for (int j = 0; j < LoadPackage.NNBndry * LoadPackage.NDim; j++)
    {
      F_e[i] += LoadPackage.M_eb[i][j] * Udd_b[j] + LoadPackage.C_eb[i][j] * Ud_b[j] + LoadPackage.K_eb[i][j] * U_b[j];
    }
  }

  // Assemble the load vector
  for (int i = 0; i < LoadPackage.NNLayer * LoadPackage.NDim; i++)
  {
    LoadPackage.UN[ND_e[i]] += F_e[i];
  }

  delete U_b;
  delete Ud_b;
  delete Udd_b;
  delete F_e;




  U_e = new double[NNLayer * NDim];
  Ud_e = new double[NNLayer * NDim];
  Udd_e = new double[NNLayer * NDim];

  F_b = new double[NNBndry * NDim];

  for (i = 0; i < NNBndry * NDim; i++)
  {
    F_b[i] = 0.0;
  }

  // Loop on the nodes on the DRM layer to find out the analytical solution (In this case only two nodes)
  for (i = 0; i < NNLayer; i++)
  {
    for (j = 0; j < NDim; j++)
    {
      x = XYZ[NoLayer_DRM[i]][j]; // Coordinate of the node
      u = v = a = 0.0;            // Initialize the values
      // Computing the analytical solution - Comment: the one-dimensional wave-motion is identical for both SV and P waves.
      if (Wave_Type == 0)
        DRM_PointValues(Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); // SV wave
      else if (Wave_Type == 1)
        DRM_PointValues(Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); // P wave

      // Filling the analytical solution vector
      U_e[i * NNBndry * NDim + j] = u;
      Ud_e[i * NNBndry * NDim + j] = v;
      Udd_e[i * NNBndry * NDim + j] = a;
    }
  }

  // Multiply the Mass, Damp, and Stiffness matrix by the vector
  for (i = 0; i < NNBndry * NDim; i++)
  {
    for (j = 0; j < NNLayer * NDim; j++)
    {
      F_b[i] += -(M_eb[j][i] * Udd_e[j] + C_eb[j][i] * Ud_e[j] + K_eb[j][i] * U_e[j]);
    }
  }

  // Assemble the load vector
  for (i = 0; i < NNBndry * NDim; i++)
  {
    UN[ND_b[i]] += F_b[i];
  }

  delete U_e;
  delete Ud_e;
  delete Udd_e;
  delete F_b;
}


/*
###################################################################################################
Purpose: This function computes ????.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 07/04/2018 - Subroutine initiated.
V1.00: 07/04/2018 - Compiled successfully.

###################################################################################################
*/

void HistorySolution(int &NJ, double &Time, double &Alpha, double &P, double &E, double &Rho, double &A, double *&U_EX, double **&XYZ)
{

  double fac;
  double x;
  double c;

  c = sqrt(E / Rho);

  fac = P * c / (Alpha * E * A);

  for (int ij = 0; ij < NJ; ij++)
  {
    x = XYZ[ij][0];
    U_EX[ij] = 0.0;

    if (Time >= (x / c))
      U_EX[ij] = U_EX[ij] + fac * (1 - cos(Alpha * (Time - x / c)));
    if (Time >= (x / c + 2 * pi / Alpha))
      U_EX[ij] = U_EX[ij] - fac * (1 - cos(Alpha * (Time - x / c - 2 * pi / Alpha)));
  }
}



/*
###################################################################################################
Purpose: This function computes 

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 07/04/2018 - Subroutine initiated.
V1.00: 07/04/2018 - Compiled successfully.

###################################################################################################
*/

void DRM_PointValues_Freq(double &amplitude, double &x, double &c, double &omega, double &u_R, double &u_I)
{

  double k; // wavenumber

  // The analytical solution is u (x,t) = u_i (exp(i k x) + exp (-i k x)) = 2 u_i cos(kx)

  k = omega / c; // wavenumber

  u_R = 2 * amplitude * cos(k * x); // The real part of the analytical solution in the frequency domain
  u_I = 0.0;                        // The imaginary part of the analytical solution in the frequency domain
}
