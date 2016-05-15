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

void der(const double&x , const vector<double>& y, 
	vector<double>& dydx)
{
	for(std::vector<double>::size_type i=0;i<y.size();++i)
		dydx[i] = -0.5*y[i]+ 0.001*y[1]  +0.2*y[(2*i)%y.size()];
}

int main()
{

	double x = 1;
	double xend = 200.;
	vector<double> y(401,1.);
	y[1] = 100;
	double h = 0.1;
		
	integrator(x,y,h,xend,eulerInt,der);


	cout << y[1] << endl << xend << endl;

	return 0;
}

