#include "NumericalOptimization.h"
#include "LIneSearchAlgorithm.h"
#include <vector>
int main()
{
	
	return 0;
}
inline double MA5171::Optimization::f(double x) {
	return std::pow((-5.0 + x), 2.0);
}
inline double MA5171::Optimization::gradient(double x) {
	return (2.0 * (-5.0 + x));
}
inline std::ofstream &MA5171::Optimization::write(std::ofstream& name,std::string fileName) {
	name.open(fileName, std::ios::app);
	return name;
}
inline double MA5171::Optimization::zoom(double alphaLo, double alphaHi) {
	double c1 = 0.0;
	double c2 = 0.0;
	double aStar = 0.0;
	int begin = 1;
	while (begin)
	{
		double aj = 0;//return interpolate result	
		if (f(aj) > f(0) + c1 * aj * gradient(0) || f(aj) >= f(alphaLo)) {
			alphaHi = aj;
		}
		else {
			if (std::fabs(gradient(aj)) <= (-c2) * gradient(0))
			{
				aStar = aj;
				break;
			}
			if ((gradient(aj)*(alphaHi-alphaLo))>=0)
			{
				alphaHi = alphaLo;

			}
			alphaLo = aj;
		}
	}
	
	return aStar;
}
double MA5171::Optimization::LineSearch(double alpha0, double alphaMax, double alpha1) {
	std::vector<double>alpha(10);
	double c1 = 0.0, c2 = 0.0;
	double astar = 0.0;
	for (size_t j = 1; j < alpha.size(); ++j)
	{
		alpha[j] = alpha[j - 1] + 0.1;
	}
	int i = 1;
	while (i) {
		if (f(alpha[i]) > f(0) + c1 * alpha[i] * gradient(0) || (f(alpha[i]) >= f(alpha[i - 1]) && i > 1)) {
			astar = zoom(alpha[i - 1], alpha[i]);
			break;
		}
		if (fabs(gradient(alpha[i])) <= -c2 * gradient(0)) {
			astar = alpha[i];
			break;
		}
		if (gradient(alpha[i])>=0.0)
		{
			astar = zoom(alpha[i], alpha[i - 1]);
			break;
		}
		++i;
	}
	return astar;
}