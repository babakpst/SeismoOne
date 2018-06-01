
#include "../include/solver_cls.h"

main_ns::solver_ns::solver_cls::solver_cls
                                (main_ns::discretization_ns::discretization_cls* aDiscretization,
                                main_ns::model_ns::model_cls* aModel):
                                DiscretizedModel(aDiscretization),
                                Model(aModel)
                                {

std::cout << " Solving the one-dimensional wave motion in the time domain ..." << std::endl;
K  = new double *[DiscretizedModel->NEqM];  // Stiffness Matrix
  for(int i=0;i<DiscretizedModel->NEqM;i++){
    K[i]=new double[DiscretizedModel->NEqM];
  }

C  = new double *[DiscretizedModel->NEqM];  // Damping matrix
  for(int i=0;i<DiscretizedModel->NEqM;i++){
    C[i]=new double[DiscretizedModel->NEqM];
  }

M  = new double *[DiscretizedModel->NEqM];  // Mass matrix
  for(int i=0;i<DiscretizedModel->NEqM;i++){
    M[i]=new double[DiscretizedModel->NEqM];
  }




K_eb  = new double *[Model->NDim * Model->NNLayer];  // 
  for(int i=0;i<(Model->NDim * Model->NNLayer);i++){
    K_eb[i]=new double[ Model->NDim * Model->NNBndry ];
  }

C_eb  = new double *[ Model->NDim * Model->NNLayer ];  // 
  for(int i=0;i<(Model->NDim * Model->NNLayer);i++){
    C_eb[i]=new double[ Model->NDim * Model->NNBndry ];
  }

M_eb  = new double *[ Model->NDim * Model->NNLayer ];  // 
  for(int i=0; i<(Model->NDim * Model->NNLayer); i++){
    M_eb[i]= new double[ Model->NDim * Model->NNBndry ];
  }

ND_b = new int  [ Model->NNBndry * Model->NDim ];
ND_e = new int  [ Model->NNLayer * Model->NDim ];

F  = new double [DiscretizedModel->NEqM] ;

// Filling the index for layered nodes
  for (int i=0;i<Model->NNLayer;i++) {
    for ( int j=0;j<Model->NDim;j++) {
      ND_e [ j * Model->NNLayer + i ] = DiscretizedModel->ID [ DiscretizedModel->NoLayer_DRM [ i ] ][j];
    }
  }

  // Filling the index for boundary nodes
  for ( int i=0;i<Model->NNBndry;i++) {
    for (int j=0;j<Model->NDim;j++) {
      ND_b [ j * Model->NNBndry + i ] = DiscretizedModel->ID [ DiscretizedModel->NoBndry_DRM[i]][j];
    }
  }

}
