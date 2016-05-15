#include "integrator.h"


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



