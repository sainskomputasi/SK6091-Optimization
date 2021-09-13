#ifndef ONE_DIMENSIONAL_IMP
#define ONE_DIMENSIONAL_IMP
#include <cmath>
#include "oneDimensional_Abstract.hpp"
double SK6091::OneDimension::goldenSec(double f(double x),double xlow, double xhigh, int maxexit, double es) {
	double xl = xlow, xu = xhigh;
	int iter = 1;
	double d = R * (xu - xl);
	double x1 = xl + d, x2 = xu - d;
	double f1 = f(x1), f2 = f(x2);
	double xopt = 0;//auxilary variable for termination 
	double fx = 0.0;//save the current value from a function 
	double ea = 0;
	if (f1 > f2)
	{
		xopt = x1;
		fx = f1;
	}
	else {
		xopt = x2;
		fx = f2;
	}
	while (iter != maxexit)
	{
		d *= R;
		if (f1 > f2)
		{
			xl = x2;
			x2 = x1;
			x1 = xl + d;
			f2 = f1;
			f1 = f(x1);
		}
		else {
			xu = x1;
			x1 = x2;
			x2 = xu - d;
			f1 = f2;
			f2 = f(x2);
		}
		++iter;
		if (f1 > f2) {
			xopt = x1;
			fx = f1;
		}
		else {
			xopt = x2;
			fx = f2;
		}
		if (xopt != 0) {
			ea = (1.0 - R) * fabs((xu - xl) / xopt) * 100;
		}
		if ((ea <= es) || (iter >= maxexit))
		{
			break;
		}
	}
	return xopt;
}
#endif