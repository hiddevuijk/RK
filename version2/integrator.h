#ifndef GUARD_integrator_h
#define GUARD_integrator_h

#include <vector>

namespace intDef {
	typedef void (*DER1)(const double&, const std::vector<double>&, std::vector<double>&);
	typedef void (*SS1)(double&,std::vector<double>&,double, DER1);
	typedef void (*DER2)(const double&, const double&, double&);
	typedef void (*SS2)(double&,double&,double, DER2);

}

void integrator(double& x, std::vector<double>& y, double h,const double& xend,
	intDef::SS1 singleStep, intDef::DER1 derivs)
{
	if(xend>x){
		while(x<xend) {
			if(xend-x<h)
				h = xend - x;
			singleStep(x,y,h,derivs);
		}
	}
	else if(xend>x) {
		while(x>xend) {
			if(xend-x>(-1*h))
				h = xend - x;
			singleStep(x,y,-1*h,derivs);
		}
	}
}


void integrator(double& x,double& y, double h,const double& xend,
	intDef::SS2 singleStep, intDef::DER2 derivs)
{
	if(xend>x){
		while(x<xend) {
			if(xend-x<h)
				h = xend - x;
			singleStep(x,y,h,derivs);
		}
	}
	else if(xend>x) {
		while(x>xend) {
			if(xend-x>(-1*h))
				h = xend - x;
			singleStep(x,y,-1*h,derivs);
		}
	}
}


#endif
