
//#include "../include/create_global_matrices_cls.h"

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

namespace main_ns
{
namespace Matrices_ns
{

class assemble_local_to_global_cls
{

  public:
    int *ND_b; // Nodal ID for DRM
    int *ND_e; // Nodal ID for DRM

    int *ND; // element constraints

    assemble_local_to_global_cls();

    virtual void assemble_local_to_global_fn(void) = 0;

    // void AssembleMassDampingStiffForceSkyline (int& NEqEl, int *& ND, double **& Ke, double **& Ce, double **& Me, double *& K_S, double *& C_S, double *& M_S, int *& JD);
};
} // namespace Matrices_ns
} // namespace main_ns
#endif
