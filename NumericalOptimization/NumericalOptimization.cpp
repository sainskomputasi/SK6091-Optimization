#include "NumericalOptimization.h"
#include "LIneSearchAlgorithm.h"
#include <Eigen/dense>
#include <vector>
#include "oneDimen/oneDimensional.hpp"
#include <chrono>
#include <cmath>
#include <math.h>
#include <corecrt_math_defines.h>
#include "specialFunctionImp.hpp"
#include "multiDimension/multiDimenImp.hpp"
double objectiveFunc(double);
template <typename T>
inline double oneD_Rastrign(const T& x);
int main()
{
	auto start = std::chrono::steady_clock::now();
	/*auto xmi = -1.5;
	SK6091::OneDimension testBrent;
	SK6091::OneDimension testGolden;
	SK6091::OneDimension testParabolic;
	double x0 = 0.0, x1 = 1.0, x2 = 4.0;
	std::cout << "maximum occurs at x \t: " <<
		testParabolic.parabolic(objectiveFunc, x0, x1, x2, 100, 0.00000000001);
	//std::cout << "minimum occurs at x \t: " << testBrent.brent(objectiveFunc, -4.0, -2.0, 0.01, 0.00000000001, &xmi, 100)<<std::endl;
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> sec = end - start;
	std::cout << "Elapsed time \t:" << sec.count()<<std::endl;
	*/
	Eigen::RowVector2d initalPoint;
	initalPoint << 0.5,0.5;
	SK6091::MultiD bfgs;
	bfgs.quasiNewton(initalPoint);
	initalPoint << 2.3, 2.3;
	bfgs.quasiNewton(initalPoint);
	initalPoint << 1.4, 2.0;
	bfgs.quasiNewton(initalPoint); 
	return 0;
}
inline double MA5171::Optimization::f(double x1,double x2) {
	return (1 + ((x1 + x2) / 4000) - (std::cos(x1) * std::cos((0.5) * x2 * std::sqrt(2.0))));
}
inline double MA5171::Optimization::fak_x1(double x1,double x2) {
	return ((x1 / 2000) + 
		(std::sin(x1) * std::cos(x2 * std::sqrt(2.0) / 2.0)));
}
inline double MA5171::Optimization::fak_x2(double x1, double x2) {
	return ((x2 / 2000) +
		(std::cos(x1) * std::sqrt(2.0) * std::sin(x2 * std::sqrt(2.0) / 2.0)) / 2.0);
}
inline std::ofstream &MA5171::Optimization::write(std::ofstream& name,std::string fileName) {
	name.open(fileName, std::ios::app);
	return name;
}
inline double MA5171::Optimization::zoom(double a_lo, double a_hi, double phi_lo, double phi_hi, double derphi_lo, double phi(double),
	double derphi(double), double phi0, double derphi0, double c1, double c2) {
	int begin = 0;
	int max = 11;
	double setFirstDel = 0.2;
	double setSectDel = 0.1; //parameter set interp ...
	double phi_rec = phi0;
	double arec = 0;
	double a, b = 0;
	double cchk,aj = 0.0; //for cub ..
	double qchk = 0.0;//for quad...
	double a_star = 0.0;
	double dAlpha = 0.0;
	while (begin!=max)
	{
		dAlpha = a_hi - a_lo;
		if (dAlpha<0)
		{
			a, b = a_hi, a_lo;
		}
		else {
			a, b = a_lo, a_hi;
		}
		if (begin > 0) {
			cchk = setFirstDel * dAlpha;
			aj = cubicInter(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, arec, phi_rec);
		}
		if ((begin==0)||(aj>b-cchk)||(aj<a+cchk))
		{
			qchk = setSectDel * dAlpha;
			aj = quadInter(a_lo, phi_lo, derphi_lo, a_hi, phi_hi);
			if ((aj==NULL)||(aj>b-qchk)||(aj<a+qchk))
			{
				aj = a_lo + 0.5 * dAlpha;//use bisec
			}
		}
		double phi_aj = phi(aj);
		//test after interpolate
		if ((phi_aj>phi0+c1*aj*derphi0)||(phi_aj>=phi_lo))
		{
			phi_rec = phi_hi;
			arec = a_hi;
			a_hi = aj;
			phi_hi = phi_aj;
		}
		else {
			double derphi_aj = derphi(aj);
			if (std::fabs(derphi_aj)<=-c2*derphi0)
			{
				a_star = aj;
				break;
			}
			if ((derphi_aj * (a_hi - a_lo))>=0) {
				phi_rec = phi_hi;
				arec = a_hi;
				a_hi = a_lo;
				phi_hi = phi_lo;
			}
			else {
				phi_rec = phi_lo;
				arec = a_lo;
			}
			a_lo = aj;
			phi_lo = phi_aj;
			derphi_lo = derphi_aj;
		}

		++begin;
	}
	return a_star;
}
/*double MA5171::Optimization::LineSearch(double alpha0, double alphaMax, double alpha1) {
	std::vector<double>alpha(10);
	double c1 = 0.0, c2 = 0.0;
	double astar = 0.0;
	for (size_t j = 1; j < alpha.size(); ++j)
	{
		alpha[j] = alpha[j - 1] + 0.1;
	}
	int i = 1;
	while (i) {
		if (f(alpha[i], alpha[i]) > f(0.0,0.0) + c1 * alpha[i] * fak_x1(0.0,0.0) || (f(alpha[i], alpha[i]) >= f(alpha[i - 1],alpha[i-1]) && i > 1)) {
			astar = zoom(alpha[i - 1], alpha[i]);
			break;
		}
		if (fabs(fak_x1(alpha[i], alpha[i])) <= -c2 * fak_x1(0.0,0.0)) {
			astar = alpha[i];
			break;
		}
		if (fak_x1(alpha[i], alpha[i])>=0.0)
		{
			astar = zoom(alpha[i], alpha[i - 1]);
			break;
		}
		++i;
	}
	return astar;
}*/

double MA5171::Optimization::cubicInter(double a, double fa, double fpa, double b, double fb, double c, double fc) {
	double C = fpa;
	double db = b - a;
	double dc = c - a;
	double denomin = std::pow((db * dc), 2.0) * (db - dc);
	//call eigen matrix template.. intended for solving the equation 
	Eigen::Matrix2d d1;
	d1 << std::pow(dc, 2.0), std::pow(-db, 2.0), 
		  std::pow(-dc, 3.0), std::pow(db, 3.0);
	Eigen::Matrix<double, 2, 1> multipl;
	multipl << (fb - fa - C * db), (fc - fa - C * dc);
	Eigen::Matrix<double, 2, 1> tempor;
	tempor << d1 * multipl;
	tempor(0, 0) /= denomin;
	tempor(1, 0) /= denomin;
	double rad = std::pow(tempor(1,0),2.0) - 3 * tempor(0, 0)*C;
	double xmin = a + (-tempor(1, 0) + std::sqrt(rad)) / (3 * tempor(0,0));
	return xmin;
}
double MA5171::Optimization::quadInter(double a, double fa, double fpa, double b, double fb) {
	double D = fa; 
	double C = fpa;
	double db = b - a * 1.0;
	double B = (fb - D - C * db) / (std::pow(db,2.0));
	double xmin = a - C / (2.0 * B); //solve for quad interp ...
	return xmin;
}
double objectiveFunc(double x) {
	return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0); 
}
double oneDrastrigin(double x) {
	return (10+std::pow(x,2.0)-10*std::cos(2*M_PI*x));
}


