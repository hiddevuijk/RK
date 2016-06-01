#ifndef GUARD_adaptRKInt_h
#define GUARD_adaptRKInt_h


#include "integrator.h"
#include <vector>

void adaptRKInt(const double& x ,double& y,const double& dy, const double& h,
		double& yerr, intDef::DER2 derivs)
{

	// perameters:
 	static const double a2=0.2,a3=0.3,a4=0.6,a5=1.,a6=0.875,
		b21=0.2,b31=3./40,b32=9./40,b41=.3,b42=-0.9,b43=1.2,
		b51=-11./54,b52=2.5,b53=-70./27,b54=35./27,
		b61=1631./55296.,b62=175./512,b63=575./13824,
		b64=44275./110592,b65=253./4096,
		c1=37./378,c2=250./621,c4=125./594,c6=512./1771;
	static const double dc1=c1-2825./27648.;
	static const double dc3=c3-18575./48384;
	static const double dc4=c4-13525./55296;
	static const double dc5=-277./14336;
	static const double dc6=c6-0.25;

	double k1,k2,k3,k4,k5,k6;

	double ytemp = y + b21*h*dy;

	derivs(x+a2*h,ytemp,k2);
	ytemp = y+h*(b31*dy+b32*k2);

	derivs(x+a3*h,ytemp,k3);
	ytemp = y+h*(b41*dy + b41*k2 + b43*k3)
	
	derivs(x+a4*h,ytemp,k4);
	ytemp = y+h*(b51*dy + b52*k2 + b53*k3 + b54*k4);

	derivs(x+a5*h,ytemp,k5);
	ytemp = y+h*(b61*dy+d62*k2+b63*k3+b64*k4+b65*k5);

	derivs(x+a6*h,ytemp,k6);
	y += h*(c1*d7+c3*k3+c4*k4+c6*k6);
	yerr = h*(dc1*dy+dc3k3+dc4*k4+dc5*k5+dc6*k6)

}





#endif
