


#ifndef SHAPEFUNCTIONS_CLS_H
#define SHAPEFUNCTIONS_CLS_H

namespace main_ns
{
  
namespace ShapeFunctions_ns
{


class ShapeFunctions_cls(){

public: 

ShapeFunctions_cls();
template<int NInt, int NNode>
void Retrieving_Gauss_Points_fn(const int& NInt, const int& NNode);

virtual void ShapeFunctions (void)  = 0;
virtual void DifferentialOfShapeFunctions (void)  = 0;
~ShapeFunctions_cls();
  
}  
}
}

#endif
