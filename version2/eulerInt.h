#ifndef GUARD_eulerInt_h
#define GUARD_eulerInt_h

#include "integrator.h"

#include <vector>

void eulerInt(double& x, std::vector<double>& y, double h, intDef::DER1 derivs);
void eulerInt(double& x, double& y, double h, intDef::DER2 derivs);

#endif
