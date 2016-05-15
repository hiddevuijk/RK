#ifndef GUARD_rkInt_h
#define GUARD_rkInt_h

#include "integrator.h"

#include <vector>

void rkInt(double& x, std::vector<double>& y, double h, intDef::DER1 derivs)
{
	std::vector<double>::size_type imax = y.size();
	std::vector<double> ym(imax,0.0);
	std::vector<double> k1(imax,0.0);
	std::vector<double> k2(imax,0.0);
	std::vector<double> k3(imax,0.0);
	std::vector<double> k4(imax,0.0);

	std::vector<double>::size_type i;
	
	derivs(x,y,k1);
	for(i=0;i<imax; ++i)
		ym[i] = y[i] + 0.5*h*k1[i];
	
	derivs(x+0.5*h,ym,k2);
	for(i=0;i<imax;++i)
		ym[i] = y[i] + 0.5*h*k2[i];

	derivs(x+0.5*h,y,k3);
	for(i=0;i<imax;++i)
		ym[i] = y[i] + h*k3[i];

	derivs(x+h,ym,k4);

	for(i=0;i<imax;++i)
		y[i] += h*(k1[i]+2*(k2[i]+k3[i])+k4[i])/6;

	x += h;
}



void rkInt(double& x,double& y, double h, intDef::DER2 derivs)
{
	double ym = 0.0;
	double k1 = 0.0;
	double k2 = 0.0;
	double k3 = 0.0;
	double k4 = 0.0;

	
	derivs(x,y,k1);
	ym = y + 0.5*h*k1;
	
	derivs(x+0.5*h,ym,k2);
	ym = y + 0.5*h*k2;

	derivs(x+0.5*h,y,k3);
	ym = y + h*k3;

	derivs(x+h,ym,k4);

	y += h*(k1+2*(k2+k3)+k4)/6;
	x += h;
}










#endif

