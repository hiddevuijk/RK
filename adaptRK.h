#ifndef GUARD_adaptRK_h
#define GUARD_adaptRK_h

#include <vector>
#include <math.h>


namespace intDef {
	typedef void (*DER1)(const double&, const std::vector<double>&, std::vector<double>&);
	typedef void (*SS1)(double&,std::vector<double>&,double, DER1);
	typedef void (*DER2)(const double&, const double&, double&);
	typedef void (*SS2)(double&,double&,double, DER2);

	typedef void (*adaptSS2)(const double&,const double&,const double&,
		const double&,double&,double&, DER2);
}


void adaptRKInt(const double& x ,const double& y,const double& dy, const double& h,
		double& ytemp, double& yerr, intDef::DER2 derivs)
{

 	static const double a2=0.2,a3=0.3,a4=0.6,a5=1.,a6=0.875,
		b21=0.2,b31=3./40,b32=9./40,b41=.3,b42=-0.9,b43=1.2,
		b51=-11./54,b52=2.5,b53=-70./27,b54=35./27,
		b61=1631./55296.,b62=175./512,b63=575./13824,
		b64=44275./110592,b65=253./4096,
		c1=37./378,c3=250./621,c4=125./594,c6=512./1771;
	static const double dc1=c1-2825./27648.;
	static const double dc3=c3-18575./48384;
	static const double dc4=c4-13525./55296;
	static const double dc5=-277./14336;
	static const double dc6=c6-0.25;

	double k2,k3,k4,k5,k6;

	ytemp = y + b21*h*dy;

	derivs(x+a2*h,ytemp,k2);
	ytemp = y+h*(b31*dy+b32*k2);

	derivs(x+a3*h,ytemp,k3);
	ytemp = y+h*(b41*dy + b42*k2 + b43*k3);
	
	derivs(x+a4*h,ytemp,k4);
	ytemp = y+h*(b51*dy + b52*k2 + b53*k3 + b54*k4);

	derivs(x+a5*h,ytemp,k5);
	ytemp = y+h*(b61*dy+b62*k2+b63*k3+b64*k4+b65*k5);

	derivs(x+a6*h,ytemp,k6);
	ytemp = y+h*(c1*dy+c3*k3+c4*k4+c6*k6);
	yerr = h*(dc1*dy+dc3*k3+dc4*k4+dc5*k5+dc6*k6);


}


void adaptIntegrator(double& x, double& y, double& xend,
	intDef::adaptSS2 singleStep, intDef::DER2 derivs)
{

	static const int maxstep = 100;
	double h = 0.01;
	static const double tiny = 1.e-30;
	double eps = 5.e-8;

	double yerr;
	double emax;
	double ytemp;
	double yscal;
	int istep = 0;

	double dy;
	while( istep<maxstep and x <= xend) {
		++istep;

		derivs(x,y,dy);
		yscal=fabs(y)+fabs(h*dy) + tiny;
		if(x+h>xend) h = xend-x;

		while(true) {	
			singleStep(x,y,dy,h,ytemp,yerr,derivs);
			emax = fabs(yerr/yscal/eps);
			if(emax<=1) break;
			h = fmax(fabs(0.9*h*pow(emax,-0.25)),0.25*fabs(h));
				
		}
		x += h;
		if(emax>1.89e-4) h=0.9*h*pow(emax,-0.2);
		else h *= 4;
		y = ytemp;
	}
}	






#endif
