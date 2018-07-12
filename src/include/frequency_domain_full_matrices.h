
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "../include/Load.h"
#include "../include/discretize_the_domain_cls.h"
#include "../include/create_global_matrices_cls.h"

#include "../include/create_full_matrices_cls.h"

#ifndef FREQUENCY_DOMAIN_FULL
#define FREQUENCY_DOMAIN_FULL

namespace main_ns
{

namespace Solver_ns
{

class frequency_domain_analysis : public main_ns::Solver_ns::Solver_cls
{


  void Reduce_the_effective_forece_in_the_freq_domain(double **&K_Eff);

  void Substitute_Freq(int &NEqM, double *&RHS, double **&K_Eff);

public:
  main_ns::discretization_ns::discretization_cls *DiscretizedModel;

  main_ns::address_ns::address_cls* Addresses;
  main_ns::model_ns::model_cls* Model;
  main_ns::discretization_ns::discretization_cls* DiscretizedModel;
  main_ns::Matrices_ns::Matrices_cls* Matrices;  

  frequency_domain_analysis(main_ns::address_ns::address_cls *, main_ns::model_ns::model_cls *,
                            main_ns::discretization_ns::discretization_cls *,
                            main_ns::Matrices_ns::Matrices_cls *);

  void Compute_the_transfer_functions_in_the_frequency_domain
  (double &alpha1, double &alpha2, int &Wave_Type, int &NEqM, double **&M, double **&C, double **&K, double **&PMat, double **&XYZ, int *&ND_e, int *&ND_b, ofstream &TransferFunc);
};

} // namespace Solver_ns
} // namespace main_ns

#endif
