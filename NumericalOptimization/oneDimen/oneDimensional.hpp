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
	std::cout << "iteration\t: " << iter << std::endl;
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
	std::cout << "Iteration \t:" << beg << std::endl;
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
		++begin;
	}
	return alpha;
}
double SK6091::OneDimension::brent(double(*f)(double), double ax, double bx, double cx, double tolerance, double* xmin, int maxIter) {
	auto begin = 1;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w,x, xm;
	double e = 0.0;
	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*f)(x);
	while (begin!=maxIter)
	{
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tolerance * std::fabs(x) + ZEPS);
		if (std::fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
			*xmin = x;
			std::cout << "total iter\t:" << begin << std::endl;
			return x;
		}
		if (std::fabs(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) {
				p = -p;
			}
			q = std::fabs(q);
			etemp = e;
			e = d;
			if (std::fabs(p) >= std::fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
				d = GOLD * (e = (x >= xm ? a - x : b - x));
			}
			else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2) {
					d = MAGN(tol1, xm - x);
				}
			}
		}
		else {
			d = GOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (std::fabs(d) >= tol1 ? x + d : x + MAGN(tol1, d));
		fu = (*f)(u);
		if (fu <= fx) {
			if (u >= x) {
				a = x;
			}
			else {
				b = x;
			}
			SHIFT(v, w, x, u);
			SHIFT(fv, fw, fx, fu);
		}
		else {
			if (u < x) {
				a = u;
			}
			else {
				b = u;
			}
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
		++begin;
	}
	*xmin = x;
	return fx;

}
#endif
