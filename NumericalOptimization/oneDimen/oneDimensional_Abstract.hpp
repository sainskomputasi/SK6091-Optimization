#ifndef ONE_DIMENSIONAL
#define ONE_DIMENSIONAL
#include <cmath>
#define R ((pow(5,0.5)-1.0)/2.0)
#define MAGN(a,b) ((b)>=0.0?std::fabs(a):-std::fabs(a))
#define GOLD 0.3819660
#define ZEPS 1.0e-10
#define SHIFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
namespace SK6091 {
	class OneDimension {
	public:
		double goldenSec(double(double), double, double, int, double);
		double parabolic(double(double), double&, double&, double&, int, double);
		double newton(double(double),double &,double,int);
		double secant(double(double), double,double, double, int);
		double brent(double(*)(double), double, double, double,double ,double *, int);

	private:
	};
}
#include "oneDimensional.hpp"
#endif