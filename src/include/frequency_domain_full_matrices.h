
#include "../include/solver_cls.h"

#ifndef FREQUENCY_DOMAIN_FULL
#define FREQUENCY_DOMAIN_FULL

namespace main_ns{

namespace Solver_ns{


class frequency_domain_analysis: public main_ns::Solver_ns::Solver_cls 
{


void LDLT_Freq ( int& NEqM, double **& K_Eff);
void Substitute_Freq ( int& NEqM, double *& RHS, double **& K_Eff);
public:


frequency_domain_analysis(main_ns::address_ns::address_cls*, main_ns::model_ns::model_cls*,
                          main_ns::discretization_ns::discretization_cls*,
                          main_ns::Matrices_ns::Matrices_cls*);

void Transfer_Full ( double & alpha1, double & alpha2, int & Wave_Type, int& NEqM, double **& M, double **& C, double **& K, double **& PMat, double **& XYZ, int *&ND_e, int *&ND_b, ofstream& TransferFunc );
};

}
}


#endif

