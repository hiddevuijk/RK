#ifndef GUARD_vectorUtil_h
#define GUARD_vectorUtil_h

#include <vector>
#include <assert.h>

void multiply(std::vector<double>& v, double x)
{
	for(int i=0;i<v.size();++i)
		v[i] *= x;
}

void devide(std::vector<double>& v, double x)
{
	assert(x != 0);
	multiply(v,1/x);
}

void add(std::vector<double>& v, double x)
{
	for(int i=0;i<v.size();++i)
		v[i] += x;
}

void subtr(std::vector<double>& v, double x)
{
	add(v,-1*x);
}

std::vector<double> add(const std::vector<double>& v1,
		const std::vector<double>& v2)
{
	int imax = v1.size() < v2.size() ? v1.size() : v2.size();
	std::vector<double> v3(imax,0.0);
	for(int i=0;i<imax;++i)
		v3[i] = v1[i] + v2[i];
	return v3;
}

std::vector<double> subtr(const std::vector<double>& v1,
		const std::vector<double>& v2)
{
	int imax = v1.size() < v2.size() ? v1.size() : v2.size();
	std::vector<double> v3(imax,0.0);
	for(int i=0;i<imax;++i)
		v3[i] = v1[i] - v2[i];
	return v3;
}




void addInplace(std::vector<double>& v1, std::vector<double>& v2)
{
	int imax = v1.size() < v2.size() ? v1.size() : v2.size();

	for(int i=0;i<imax;++i)
		v1[i] += v2[i];
}

void subtrInplace(std::vector<double>& v1, std::vector<double>& v2)
{
	int imax = v1.size() < v2.size() ? v1.size() : v2.size();
	for(int i=0;i<imax;++i)
		v1[i] -= v2[i];
}






#endif
