#ifndef GUARD_eulerInt_h
#define GUARD_eulerInt_h

#include "integrator.h"

#include <vector>

void eulerInt(double& x, std::vector<double>& y, double h, intDef::DER1 derivs)
{
	std::vector<double> dydx(y.size(),0.0);
	derivs(x,y,dydx);
	for(std::vector<double>::size_type i=0;i<y.size();++i)
		y[i] += h*dydx[i];
	x += h;
}


void eulerInt(double& x, double& y, double h, intDef::DER2 derivs)
{
	double dydx = 0.0;
	derivs(x,y,dydx);
	y += h*dydx;
	x += h;
}



#endif
