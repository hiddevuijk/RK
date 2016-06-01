#ifndef GUARD_adaptRKInt_h
#define GUARD_adaptRKInt_h

#include <vector>

void adaptRKInt(const double& y, const double& dy,
	const double& x, const double& h, double& yout,
	double& yerr)
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

	double ytemp = y + b21*h*dy;
	derivs(x+a2*h,ytemp,k2);

	ytemp = y+h*(b31*dy+b32*k2);
	derivs(






#endif
