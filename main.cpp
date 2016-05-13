/*

main program for testing the Runge-Kutta integration scheme

*/


#include "integrator.h"
#include "eulerInt.h"


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
		dydx[i] = 1.*i + y[i]*x;
}

int main()
{

	double x = 1;
	double xend = 10.;
	vector<double> y(2,1.);
	double h = 0.1;
	
	integrator(x,y,h,xend,eulerInt,der);


	cout << x << endl << xend << endl;

	return 0;
}

