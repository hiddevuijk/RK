#ifndef GUARD_integrator_h
#define GUARD_integrator_h

#include <vector>

typedef void (*DER)(const double&, const std::vector<double>&, std::vector<double>&);
typedef void (*SS)(double&,std::vector<double>&,double, DER);


void integrator(double& x, std::vector<double>& y, double h,const double& xend,
	SS singleStep, DER derivs)
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
