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
double SK6091::OneDimension::parabolic(double f(double), double &x0, double &x1, double &x2, int maxiter, double tolerance) {
	double x3 = 0.0;
	int beg = 0;
	std::cout << "test 1:\t: "<< f(1)<<std::endl;

	while (std::fabs(x2-x3)>tolerance) {
		x3 = (f(x0) * (std::pow(x1, 2.0) - std::pow(x2, 2.0)) + 
			f(x1) * (std::pow(x2, 2.0) - std::pow(x0, 2.0)) + f(x2) * (std::pow(x0, 2.0) - std::pow(x1, 2.0))) / (2.0*f(x0)*(x1-x2) +
			2.0 * f(x1) * (x2 - x0) + 2.0 * f(x2) * (x0 - x1));
		if ((f(x1) < f(x3)) && (x3>x1)) {
			x0 = x1;
			x1 = x3;
		}
		else {
			x2 = x1;
			x1 = x3;
		}
		++beg;

	}
	return x3;
}
double SK6091::OneDimension::newton(double f(double x), double &guess, double tolerance, int maxIter) {
	auto x0 = guess;
	double x1 = 0.0;
	int begin = 1;
	auto delx = 0.01;
	auto xi = 0.0;
	auto secDerivative = 0.0,firstDerivative=0.0;
	while (begin!=maxIter)
	{
		firstDerivative = (f(guess + delx) - f(guess - delx)) / (2.0 * delx);
		secDerivative = (f(guess + delx) + f(guess - delx) - 2 * f(guess)) / (std::pow(delx, 2.0));
		xi = guess;
		guess = xi - firstDerivative / secDerivative;
		if (std::fabs(guess - xi) < tolerance) {
			break;
		}
		++begin;
	}
	return guess;
}
double SK6091::OneDimension::secant(double f(double), double lo, double hi, double tolerance, int maxIter) {
	auto sign = 0;
	auto delx = 0.01;
	int begin = 1;
	auto alpha = 0.0;
	auto dervA = 0.0;
	auto dervB = 0.0;
	auto derAlpha = 0.0;
	while (begin!=maxIter)
	{
		alpha = (lo + hi) / 2.0;
		dervA = (f(lo + delx) - f(lo - delx)) / (2.0 * delx);
		derAlpha = (f(alpha + delx) - f(alpha - delx)) / (2.0 * delx);
		dervB = (f(hi + delx) - f(hi - delx)) / (2.0 * delx);
		if ((dervA*derAlpha)<0) {
			hi = alpha;
			sign = 1;
		}
		else {
			lo = alpha;
		}
		if (sign==1)
		{
			alpha = hi - (dervB / ((dervB - dervA) / (hi - lo)));
			derAlpha = (f(alpha + delx) - f(alpha - delx)) / (2.0 * delx);
			if (derAlpha>0)
			{
				hi = alpha;
			}
			else {
				lo = alpha;
			}
			if (std::fabs(derAlpha) < tolerance)
			{
				break;
			}
		}
	}
	return alpha;
}
#endif
