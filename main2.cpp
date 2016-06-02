/*

main program for testing the Runge-Kutta integration scheme

*/


#include "integrator.h"
#include "eulerInt.h"
#include "rkInt.h"
#include "adaptRKInt.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace::std;

void der(const double&x , const double& y, double& dydx)
{
	dydx = -0.1*y + y/x;

}

int main()
{

	double x = 1;
	double xend = 20.;
	double y = 1.;
	adaptIntegrator(x,y,xend,adaptRKInt,der);


	cout << y << endl << xend << endl;

	return 0;
}

