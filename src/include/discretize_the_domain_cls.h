
// include files
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "../include/reading_the_model_cls.h"

#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

namespace main_ns
{
namespace discretization_ns
{

class discretization_cls
{

  // members
private:
  const main_ns::model_ns::model_cls *model;

public:
  int NEqM; // number of Equations (Modified i.e. after applying essential boundary conditions)

  int *MTel;  // material type of each Element
  int **INod; // node connectivity
  int **ID;   // identity

  int *Nodal_History; // a vector that holds the node no. to record the history of displc
  int *NoBndry_DRM;   // a vector that holds the node numbers on the DRM boundary
  int *NoLayer_DRM;   // a vector that holds the node numbers on the DRM layer

  double **XYZ; // node coordinates

  // methods
public:
  explicit discretization_cls(const main_ns::model_ns::model_cls *);
  void Discretization();

}; // class discretization_cls
} // namespace discretization_ns
} // namespace main_ns

#endif
