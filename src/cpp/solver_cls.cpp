
#include "../include/solver_cls.h"

main_ns::Solver_ns::Solver_cls::Solver_cls(
                            main_ns::discretization_ns::discretization_cls* aDiscretization, 
                            main_ns::model_ns::model_cls* aModel,
                            main_ns::Matrices_ns::Matrices_cls* aMatrices)
                            : DiscretizedModel(aDiscretization), Model(aModel), Matrices(aMatrices)
{
  
}