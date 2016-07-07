/*

main program for testing the Runge-Kutta integration scheme

*/


#include "integrator.h"
#include "eulerInt.h"
#include "rkInt.h"

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
	double h = 1.e-5;	
	integrator(x,y,h,xend,rkInt,der);


	cout << y << endl << xend << endl;

	return 0;
}





