#ifndef GUARD_integrate_h
#define GUARD_integrate_h

#include <vector>
#include <math.h>
#include <iostream>

namespace intDef {
	typedef void (*DER1)(const double&, const std::vector<double>&, std::vector<double>&);
	typedef void (*SS1)(double&,std::vector<double>&,double, DER1);
	typedef void (*DER2)(const double&, const double&, double&);
	typedef void (*SS2)(double&,double&,double, DER2);

	typedef void (*adaptSS2)(const double&,const double&,const double&,
		const double&,double&,double&, DER2);
	typedef void (*adaptSS1 )(const double&,const std::vector<double>&,
		const std::vector<double>&, const double&,std::vector<double>&, 
		double&,intDef::DER1);

	double maxVec(const std::vector<double>&);
}


void adaptRKstep(const double& x ,const double& y,const double& dy, const double& h,
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

void adaptRKstep(const double& x ,const std::vector<double>& y,
		const std::vector<double>& dy,const double& h,
		std::vector<double>& ytemp, double& yerr,intDef::DER1 derivs)
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

	std::vector<double>::size_type imax = y.size();
	std::vector<double> yerrVec(imax,0.0);
	std::vector<double> k2(imax,0.0);
	std::vector<double> k3(imax,0.0);
	std::vector<double> k4(imax,0.0);
	std::vector<double> k5(imax,0.0);
	std::vector<double> k6(imax,0.0);

	std::vector<double>::size_type i;

	for(i=0;i<imax;++i)
		ytemp[i] = y[i] + b21*h*dy[i];

	derivs(x+a2*h,ytemp,k2);
	for(i=0;i<imax;++i)
		ytemp[i] = y[i]+h*(b31*dy[i]+b32*k2[i]);

	derivs(x+a3*h,ytemp,k3);
	for(i=0;i<imax;++i)
		ytemp[i] = y[i]+h*(b41*dy[i]+b42*k2[i]+b43*k3[i]);
	
	derivs(x+a4*h,ytemp,k4);
	for(i=0;i<imax;++i)
		ytemp[i] = y[i]+h*(b51*dy[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);

	derivs(x+a5*h,ytemp,k5);
	for(i=0;i<imax;++i)
		ytemp[i] = y[i]+h*(b61*dy[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);

	derivs(x+a6*h,ytemp,k6);
	for(i=0;i<imax;++i)
		ytemp[i] = y[i]+h*(c1*dy[i]+c3*k3[i]+c4*k4[i]+c6*k6[i]);

	yerr = dc1*dy[0]+dc3*k3[0]+dc4*k4[0]+dc5*k5[0]+dc6*k6[0];
	for(i=1;i<imax;++i) {
		double erri = dc1*dy[i]+dc3*k3[i]+dc4*k4[i]+dc5*k5[i]+dc6*k6[i];
		yerr = yerr > erri ? yerr : erri;
	}	
	yerr *= h;



}

/*
 * euler step
 */

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



/*
	Fourth-order Runge-Kutta integration, with multiple
	dependent variables.

	Independent variable:	x
	Dependent variables:	y[i]
	Stepsize:	h
	Derivatives are defined in derivs; see integrator.h.
*/

void RKstep(double& x, std::vector<double>& y, double h, intDef::DER1 derivs)
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

/*
	Fourth-order Runge-Kutta integration, with a single dependent variable.

	Independent variable:	x
	Single dependent variables:	y
	Stepsize:	h
	Derivatives are defined in derivs; see integrator.h.
*/



void RKstep(double& x,double& y, double h, intDef::DER2 derivs)
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


/* *************************************************
 *
 *  integration methods
 *
 *  ************************************************/

void adaptIntegrator(double& x, double& y, double& xend,
	intDef::adaptSS2 singleStep, intDef::DER2 derivs,
	double h, const double& eps,
	const unsigned long int& maxstep)
{

	static const double tiny = 1.e-30;

	double yerr;
	double emax;
	double ytemp;
	double yscal;
	unsigned long int istep = 0;

	double dy;
	while( istep<maxstep and x < xend) {
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

void adaptIntegrator(double& x, std::vector<double>& y, double& xend,
	intDef::adaptSS1 singleStep, intDef::DER1 derivs,
	double h, const double& eps,
	const unsigned long int& maxstep)
{
	static const double tiny = 1.e-30;

	double yerr;
	double emax;
	double yscal;
	unsigned long int istep = 0;

	std::vector<double>::size_type i;
	std::vector<double>::size_type imax = y.size();
	std::vector<double> dy(imax);
	std::vector<double> ytemp(imax,0.0);
	while( istep<maxstep and x <= xend) {
		++istep;
		derivs(x,y,dy);
		yscal = intDef::maxVec(y)+h*intDef::maxVec(dy)+tiny;
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

void Integrator(double xi, const double yi, double )

/*
 *  Auxiliary functions
 */

double intDef::maxVec(const std::vector<double>& v)
{
	double max = v[0];
	std::vector<double>::const_iterator it = v.begin();
	while(it != v.end()) {
		if (*it > max) max = *it;
		++it;
	}
	return max;
}


#endif
