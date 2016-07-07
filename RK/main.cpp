/*

main program for testing the Runge-Kutta integration scheme

*/


#include "adaptRK.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace::std;

void der(const double&x , const vector<double>& y, vector<double>& dydx)
{
	for(int i=0;i<y.size();++i)
		dydx[i] = -0.1*y[i] + y[i]/x ;

}
void der0(const double&x , const double& y, double& dydx)
{

	dydx = -0.1*y+ y/x;

}


int main()
{
	double eps = 1e-15;
	double h0 = 0.001;
	double x = 1.;
	double xend = 20.;
	vector<double> y(2,1.);
	adaptIntegrator(x,y,xend,adaptRKInt,der,h0,eps,1e7);

	x = 1.;
	double yy = 1.;
	adaptIntegrator(x,yy,xend,adaptRKInt,der0,h0,eps,1e7);
	cout << y[0] << '\t' << y[1]  << endl;
	cout << yy << endl;
	cout << endl << y[0] - yy << endl;
	return 0;
}

