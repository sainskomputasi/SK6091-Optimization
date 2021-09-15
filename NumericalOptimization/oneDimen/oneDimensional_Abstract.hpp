#ifndef ONE_DIMENSIONAL
#define ONE_DIMENSIONAL
#include <cmath>
#define R ((pow(5,0.5)-1.0)/2.0)
namespace SK6091 {
	class OneDimension {
	public:
		double goldenSec(double(double), double, double, int, double);
		double parabolic(double(double), double&, double&, double&, int, double);
		double newton(double(double),double &,double,int);
	private:
	};
}
#include "oneDimensional.hpp"
#endif