#ifndef GUARD_rkInt_h
#define GUARD_rkInt_h

#include "integrator.h"

#include <vector>

void rkInt(double& x, std::vector<double>& y, double h, intDef::DER1 derivs);
void rkInt(double& x,double& y, double h, intDef::DER2 derivs);

#endif

