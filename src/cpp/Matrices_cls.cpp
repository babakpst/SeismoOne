
#include "../include/Matrices_cls.h"

main_ns::Matrices_ns::Matrices_cls::Matrices_cls
                                (main_ns::discretization_ns::discretization_cls* aDiscretization,
                                main_ns::model_ns::model_cls* aModel):
                                DiscretizedModel(aDiscretization),
                                Model(aModel)
                                {}
