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
#include <random>
void solveNonLinear(Eigen::Matrix<double,1,4>,Eigen::Matrix<double,1,11>,Eigen::Matrix<double,1,11>);
double objectiveFunc(double);
template <typename T>
inline double oneD_Rastrign(const T& x);
int main()
{
	//auto start = std::chrono::steady_clock::now();
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
	
	/*Eigen::RowVector2d initalPoint;
	initalPoint << 0.5,-0.5;
	SK6091::MultiD bfgs;
	bfgs.stepestDes(initalPoint);
	initalPoint << 0.1, 0.1;
	bfgs.stepestDes(initalPoint);
	initalPoint << -0.1, -0.1;
	bfgs.stepestDes(initalPoint);
	initalPoint << 0.6, 0.9;
	bfgs.stepestDes(initalPoint);
	initalPoint << 4.9, 4.9;
	bfgs.stepestDes(initalPoint); 
	*/
	/*auto start = std::chrono::steady_clock::now();
	Eigen::RowVector2d initalPoint;
	initalPoint << 0.5,-0.5;
	SK6091::functionTest::hessian(initalPoint);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> sec = end - start;
	std::cout << "Elapsed time \t:" << sec.count() << std::endl; */
	//implementation newton system  
	/*auto begin = 1;
	Eigen::RowVector3d X,Y,ftemp;
	X << 0.1, 0.1, -0.1;
	Eigen::Matrix3d A;
	auto tol = 0.0000001;
	A << 0, 0, 0, 0, 0, 0, 0, 0, 0;
	std::cout << "x(0)\t x(1) \t x(2)" << std::endl;
	while (begin!=100)
	{
		for (size_t i = 0; i <= 2; i++)
		{
			for (size_t j = 0; j <= 2; j++) {
				A(i , j ) = SK6091::functionTest::P(i+1, j+1, X); //3*3 for jacobian
			}
		}
		for (size_t i = 0; i <= 2; i++)
		{
			ftemp(i) = -SK6091::functionTest::F(i+1, X);// 3*1 for -f
		} 
		
		Y = (A.inverse() * ftemp.transpose()).transpose();
		X = X + Y;
		std::cout << X[0] << "\t" << X[1] << "\t" <<X[2] << std::endl;
		if (Y.norm() < tol)
		{
			std::cout << "Coverge to .....\t: " << std::endl;
			std::cout << X << std::endl;
			break;
		}
		++begin;
	}
	*/
	//test for broyden method.... 
	/*auto begin = 1;
	Eigen::RowVector3d u,X,z,Y, s, v,vv,w;
	X << 0.1, 0.1, -0.1;
	auto p = 0.0;
	Eigen::Matrix3d A,A1;
	auto tol = 0.0000001;
	A << 0, 0, 0, 0, 0, 0, 0, 0, 0;
	std::cout << "x(0)\t x(1) \t x(2)" << std::endl;
	for (size_t i = 0; i <= 2; i++)
	{
		for (size_t j = 0; j <= 2; j++) {
			A(i, j) = SK6091::functionTest::P(i + 1, j + 1, X); //3*3 for jacobian
		}
	}
	for (size_t i = 0; i <= 2; i++)
	{
		v(i) = SK6091::functionTest::F(i + 1, X);// 3*1 for -f
	}
	A1 = A.inverse();
	s = (-A1 * v.transpose()).transpose();
	X += s;
	auto k = 2;
	while (k<100)
	{
		w = v;

		for (size_t i = 0; i <= 2; i++)
		{
			v(i) = SK6091::functionTest::F(i + 1, X);// 3*1 for -f
		}
		Y = v - w;

		z = (-A1 * Y.transpose()).transpose();
		p = -s * z.transpose();
		u = s * A1;
		A1 = A1 +  (1.0/p)*(s + z).transpose() * u;
		s = (-A1 * v.transpose()).transpose();
		X = X + s;
		if (s.norm()<0.00001)
		{
			std::cout << "coverge to .....  with total iteration ... " << k - 1 << std::endl;;
			std::cout << X << std::endl;
			break;
		}
		++k;
	}
	*/
	//test for steepest descent ...
	/*auto k = 1;
	const int max = 100;
	auto z0 = 0.0;
	Eigen::RowVector3d g1,X;
	X << 0.1, 0.1, -0.1;
	Eigen::RowVector4d G;
	auto g0 = 0.0, g3 = 0.0, g2 = 0.0;;
	Eigen::RowVector3d z;
	auto alpha0=0.0,alpa1 = 0.0, alpha3 = 1.0,alpha2 =0.0;
	auto h1 = 0.0, h2 = 0.0, h3 = 0.0;
	auto a0 = 0.0;
	while (k<max)
	{
		G[0] = SK6091::functionTest::sumG(X, 3);
		for (size_t i = 0; i <3; i++)
		{
			z[i] = SK6091::functionTest::sumDerph(i + 1, X);
		}
		z0 = z.norm();
		if (z0==0)
		{
			std::cout << "Zero gradient ..... \n" << X << std::endl;
		}
		//assign here .... 
		z = z / z0;
		g3 = SK6091::functionTest::sumG(X-alpha3*z, 3);

		while (g3>=G[0])
		{
			alpha3 = alpha3 / 2.0;
			if (alpha3< 1.0E-20)
			{
				std::cout << "No likely improvement ...." << std::endl;	
				std::cout << X << std::endl;
				break;
			}
			else {
				g3 = SK6091::functionTest::sumG(X - alpha3 * z, 3);
			}
		}
		//prepare for interpolation.... 
		alpha2 = alpha3 / 2.0;
		g2 = SK6091::functionTest::sumG(X - alpha2 * z, 3);
		h1 = (g2 - G[0]) / alpha2;
		h2 = (g3 - g2) / (alpha3 - alpha2);
		h3 = (h2 - h1) / alpha3;
		alpha0 = 0.5 * (alpha2 - h1 / h3);
		g0 = SK6091::functionTest::sumG(X - alpha0 * z,3);
		a0 = alpha0;
		
			if (std::fabs(G[0]) < std::fabs(g0)) {
				a0 = alpha0;
				g0 = G[0];
			}
			if (std::fabs(g2) < std::fabs(g0)) {
				a0 = alpha0;
				g0 = G[1];
			}
			if (std::fabs(g3) < std::fabs(g0)) {
				a0 = alpha0;
				g0 = G[2];
			}
		
		X = X - a0 * z;
		if ((g0< 1.0E-20) || (std::fabs(g0 - G[0])< 1.0E-20))
		{
			std::cout <<"coverge , with total iteration \t:" << k << std::endl;
			std::cout << "value\n" << X << std::endl;
			break;
		}
		++k;
	}
	*/
	std::cout << "test for non linear least square" << std::endl;
	auto x1 = 2.54107 , x2 = 0.259502;
	Eigen::Matrix<double, 1, 11> t_i;
	t_i << 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0;
	Eigen::Matrix<double, 1, 11> y_i;
	/*std::cout << "test coverge to solution \t: [";
	for (size_t i = 0; i < 5; i++)
	{
		std::cout << x1 * std::exp(x2 * t_i(0, i)) <<"| ";
	}
	std::cout << "]\n"; */
	//y_i << 3.2939, 4.2699, 7.1749, 9.3008, 20.259;
	y_i << 0.0, 3.55, 3.82, 2.98, 2.32, 1.48, 1.02, 0.81, 0.41, 0.42, 0.15;
	//13.6189 - 0.208505 - 13.6218 - 0.462183
	//-13.6219 -0.462181    13.619 -0.208506
	Eigen::Matrix<double, 1, 4> X;
	X << -12.6218, -0.362183, 13.6189, -0.208505;
	solveNonLinear(X, y_i, t_i);
	std::cout<<"Result: "<<SK6091::functionTest::gaussNewton(X, y_i, t_i);

	X << 13.6189, - 0.208505, - 13.6218 ,- 0.462183;
	solveNonLinear(X, y_i, t_i);
	std::cout << "Prove solution .... " << std::endl;
	X << -13.6218, -0.462183, 13.6189, -0.208505;

	for (size_t i = 0; i < 11; i++)
	{
		std::cout << X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i)) <<"| ";
	}
	std::cout << "1.experimant model ..." << std::endl;
	X << -13.5218, -0.462183, 13.6189, -0.208505;
	for (size_t i = 0; i < 11; i++)
	{
		std::cout << X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i)) << "| ";
	}
	std::cout << "2,experimant model ..." << std::endl;

	X << 5.68603 ,- 0.124843 ,- 5.6543 ,- 0.842367;
	for (size_t i = 0; i < 11; i++)
	{
		std::cout << X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i)) << "| ";
	}
	X << 14.6189, -0.108505, -13.6218, -0.262183;
	X << -12.6218, -0.362183, 13.6189, -0.208505;
	std::cout << "\ntest for levenberg method \t:\n";
	//std::cout << SK6091::functionTest::Levenberg(X, y_i, t_i) << std::endl;
	X << -6.10177, -0.902389, 6.15153, -0.135121;
	X << 14.6189, -0.108505, -13.6218, -0.262183;

	//std::cout << SK6091::functionTest::Levenberg(X, y_i, t_i) << std::endl;
	X << 5.68603, -0.124843, -5.6543, -0.842367;
	std::cout << "3. experimant model ..." << std::endl;

	for (size_t i = 0; i < 11; i++)
	{
		std::cout << X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i)) << "| ";
	}
	//test for dogleg method  .... 
	//std::cout << "Test for dogled method ...." << std::endl;
	//std::cout << SK6091::functionTest::dogLeg(X, y_i, t_i) << std::endl;
	/*std::cout << "test homotopy \t:" << std::endl;
	Eigen::RowVector3d initialCondition;
	initialCondition << 0.0, 0.0, 0.0;
	SK6091::functionTest::homotopy(initialCondition);
	std::cout << "\ntest for n =2" << std::endl;
	Eigen::RowVector2d initialCon;
	initialCon << 0.738688 ,0.673325;
	SK6091::functionTest::homotopyt(initialCon);
	std::cout << "debug equation 1\t: " <<std::sin(0.739085 ) - 0.673612 <<std::endl;
	std::cout <<"debug equation 2\t: "<< std::pow(0.739085, 2.0) + std::pow(0.673612, 2.0) - 1.0 << std::endl;
	*/
	/*std::cout << "Polak ribiere implementation " << std::endl;
	Eigen::RowVector2d initialGuessPolak;
	initialGuessPolak << 0.8, 0.9;
	SK6091::functionTest::polakRibere(initialGuessPolak);
	initialGuessPolak << 1.0, 1.0;
	SK6091::functionTest::polakRibere(initialGuessPolak);
	initialGuessPolak << 0.1, 0.1;
	SK6091::functionTest::polakRibere(initialGuessPolak);
	*/
	//test for levenberg
	/*Eigen::RowVector2d initialGuessPolak;
	initialGuessPolak << 1.0, 0.5;
	SK6091::functionTest::levenvergMarquardt(initialGuessPolak); */
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
void solveNonLinear(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i) {
	auto begin = 1, end = 100;
	Eigen::Matrix<double, 1, 4> p;
	p << 0, 0, 0, 0;
	while (begin != end)
	{
		p = SK6091::functionTest::hes(X, y_i, t_i).inverse() * (-SK6091::functionTest::fDerp(X, y_i, t_i));
		if (p.norm() < 1.e-15)
		{
			std::cout << "coverge to :\n";
			std::cout << X << std::endl;
			break;
		}
		X = X + p;
		++begin;
	}
}


