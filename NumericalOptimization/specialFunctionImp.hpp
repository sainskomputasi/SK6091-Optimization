#ifndef SPECIAL_FUNCTION_IMP
#define SPECIAL_FUNCTION_IMP
#include <cmath>
#include <math.h>
#include <eigen/dense>
#include "specialFunction.hpp"
inline double SK6091::functionTest::Griewank(Eigen::RowVector2d temp){
	return (1 + (1 / 4000) * std::pow(temp[0], 2.0) + (1.0 / 4000.0) * std::pow(temp[1], 2.0) - std::cos(temp[0]) * std::cos((1.0 / 2.0) * temp[1] * std::sqrt(2.0)));
}
inline  double SK6091::functionTest::Rosenbrock(Eigen::RowVector2d temp) {
	return (100 * (std::pow(temp[1] - std::pow(temp[0], 2.0), 2.0)) + std::pow(1 - temp[0], 2.0));
}
inline double SK6091::functionTest::Ackley(Eigen::RowVector2d temp) {
	double a = 20.0, b = 0.2, c = 2.0 * M_PI;
	auto begin = 0, end = 2;
	auto s1 = 0.0, s2 = 0.0;
	while (begin!=end)
	{
		s1 = s1 + std::pow(temp[begin], 2.0);
		s2 = s2 + std::cos(c * temp[begin]);
		++begin;
	}
	return (-a * std::exp(-b * std::sqrt(1.0 / end * s1)) - std::exp(1.0 / end * s2) + a + std::exp(1.0));
}
inline double SK6091::functionTest::Booth_s(Eigen::RowVector2d temp) {
	auto s = 0.0;
	s = std::pow((temp[0] + 2 * temp[1] - 7.0), 2.0) + std::pow((2 * temp[0] + temp[1] - 5.0),2.0);
	return s;
}
inline double SK6091::functionTest::Michalewicz(Eigen::RowVector2d temp) {
	auto begin = 0, end = 2;
	auto s = 0.0;
	auto m = 10.0;
	while (begin!=end)
	{
		s = s + std::sin(temp[begin]) *
			(std::pow(
				std::sin( ((begin+1) * std::pow(temp[begin], 2.0) / M_PI) ), (2.0 * m) )
				);
		++begin;
	}
	return (-s);
}
inline double SK6091::functionTest::Spring_sys(Eigen::RowVector2d temp) {
	return (100 * std::pow(std::sqrt(std::pow(temp[0], 2.0) + std::pow((temp[1] + 1.0), 2.0)) - 1.0, 2.0) + 90 * std::pow(std::sqrt(std::pow(temp[0], 2.0) + std::pow((temp[1] - 1.0), 2.0)) - 1.0, 2.0) - (20 * temp[0] + 40 * temp[1]));
}
inline double SK6091::functionTest::Rastrigin(Eigen::RowVector2d temp) {
	auto begin = 0, end = 2;
	auto s = 0.0;
	while (begin!=end)
	{
		s = s + (std::pow(temp[begin], 2.0) - 10 * std::cos(2 * M_PI * temp[begin]) );
		++begin;
	}
	return (10*2+s);
}
inline Eigen::RowVector2d SK6091::functionTest::grad(Eigen::RowVector2d vecName) {
	Eigen::RowVector2d temp;
	temp << 0, 0;
	auto xvec = vecName;
	auto delx = 1e-3;
	auto xvec1 = vecName;
	for (size_t i = 0; i != 2; ++i)
	{
		xvec = vecName;
		xvec1 = vecName;
		xvec[i] = vecName[i] + delx;
		xvec1[i] = vecName[i] - delx;
		temp[i] = (Griewank(xvec) - Griewank(xvec1)) / (2 * delx);
	}
	return  temp;
}
Eigen::Matrix2d SK6091::functionTest::hessian(Eigen::RowVector2d vecName) {
	auto del = 1e-3;
	Eigen::Matrix2d temporary(2, 2);
	auto temp = vecName;
	temp << 0, 0;
	double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
	temporary << 0, 0, 0, 0;
	std::ofstream write;
	write.open("hessian.csv", std::ios::app);
	for (size_t i = 0; i < 2; ++i)
	{
		for (size_t j = 0; j < 2; ++j) {
			if (i == j) {
				temp = vecName;
				temp[i] = vecName[i] + del;
				term1 = Rastrigin(temp);
				temp[i] = vecName[i] - del;
				term2 = Rastrigin(temp);
				term3 = Rastrigin(vecName);
				temporary(i, j) = (term1 - 2.0 * term3 + term2) / (std::pow(del, 2.0));
			}
			else {
				temp = vecName;
				temp[i] = vecName[i] + del;
				temp[j] = vecName[j] + del;
				term1 = Rastrigin(temp);

				temp = vecName;
				temp[i] = vecName[i] + del;
				temp[j] = vecName[j] - del;
				term2 = Spring_sys(temp);

				temp = vecName;
				temp[i] = vecName[i] - del;
				temp[j] = vecName[j] + del;
				term3 = Rastrigin(temp);


				temp = vecName;
				temp[i] = vecName[i] - del;
				temp[j] = vecName[j] - del;
				term4 = Rastrigin(temp);
				temporary(i, j) = (term1 - term2 - term3 + term4) / (4.0 * std::pow(del, 2.0));
			}
		}
	}
	write << temporary(0, 0) << "|" << temporary(0, 1) << "|"<<temporary(1, 0) << "|"<<temporary(1, 1) << std::endl;
	write.close();
	return temporary;

}
Eigen::RowVector2d SK6091::functionTest::goldFunc(Eigen::RowVector2d x, Eigen::RowVector2d search) {
	double a = -5.0, b = 5.0;
	auto tau = 0.381967;;
	double alpha1 = a * (1 - tau) + b * tau;
	double alpha2 = a * tau + b * (1 - tau);
	double falpha1 =	Griewank(x + alpha1 * search);
	double falpha2 = Griewank(x + alpha2 * search);
	auto begin = 1;
	double epsilon = 1e-5;
	Eigen::RowVector2d temp;
	temp << 0, 0;
	while (begin != 1000)
	{
		if (falpha1 > falpha2) {
			a = alpha1;
			alpha1 = alpha2;
			falpha1 = falpha2;
			alpha2 = tau * a + (1 - tau) * b;
			falpha2 = Griewank(x + alpha2 * search);
		}
		else {
			b = alpha2;
			alpha2 = alpha1;
			falpha2 = falpha1;
			alpha1 = tau * b + (1 - tau) * a;
			falpha1 = Griewank(x + alpha1 * search);
		}
		if (std::fabs(Griewank(x + alpha1 * search) - Griewank(x + alpha2 * search)) < epsilon) {
			temp[0] = alpha1;
			temp[1] = falpha1;
			return temp;
			break;
		}
		++begin;
	}
	return temp;
}
inline double SK6091::functionTest::F(int I, Eigen::RowVector3d X) {
	double f = 0.0;

	switch (I) {
	case 1:
		f = 3 * X[0] - cos(X[1] * X[2]) - 0.5;
		break;
	case 2:
		f = X[0] * X[0] - 81 * (X[1] + 0.1) * (X[1] + 0.1) + sin(X[2]) + 1.06;
		break;
	case 3:
		f = exp(-X[0] * X[1]) + 20 * X[2] + (10 * M_PI - 3) / 3;
		break;
	}
	return f;

}
inline  double SK6091::functionTest::P(int I, int J, Eigen::RowVector3d X) {
	double p = 0.0;

	switch (I) {
	case 1:
		switch (J) {
		case 1:
			p = 3;
			break;
		case 2:
			p = X[2] * sin(X[1] * X[2]);
			break;
		case 3:
			p = X[1] * sin(X[1] * X[2]);
			break;
		}
		break;
	case 2:
		switch (J) {
		case 1:
			p = 2 * X[0];
			break;
		case 2:
			p = -162 * (X[1] + 0.1);
			break;
		case 3:
			p = cos(X[2]);
			break;
		}
		break;
	case 3:
		switch (J) {
		case 1:
			p = -X[1] * exp(-X[0] * X[1]);
			break;
		case 2:
			p = -X[0] * exp(-X[0] * X[1]);
			break;
		case 3:
			p = 20;
			break;
		}
		break;
	}
	return p;

}
inline  double SK6091::functionTest::CF(int i, Eigen::RowVector3d X) {
	double cf = 0.0;
	switch (i)
	{
		case 1:
			cf = 3 * X[0] - cos(X[1] * X[2]) - 0.5;
			break;
		case 2:
			cf = X[0] * X[0] - 81 * (X[1] + 0.1) * (X[1] + 0.1) + sin(X[2]) + 1.06;
			break;
		case 3:
			cf = exp(-X[0] * X[1]) + 20 * X[2] + (10 * M_PI - 3) / 3;
			break;
	}
	return cf;
}
inline  double  SK6091::functionTest::sumDerph(int i, Eigen::RowVector3d X) {
	double p = 0.0; 
	switch (i)
	{
	case 1:
		p = 2*(3)*CF(1,X) + 
			2*(2 * X[0]) * CF(2,X) + 
			2* (-X[1] * exp(-X[0] * X[1]))* CF(3,X);
		break;
	case 2:
		p = 2 * (X[2] * sin(X[1] * X[2])) * CF(1, X) +
			2 * (-162 * (X[1] + 0.1)) * CF(2, X) +
			2 * (-X[0] * exp(-X[0] * X[1])) *CF(3,X);
		break;
	case 3:
		p = 2 * (X[1] * sin(X[1] * X[2])) * CF(1, X)
			+ 2 * (cos(X[2])) * CF(2, X) +
			2 * 20 * CF(3, X);
		break; 
	}
	return p;
}
inline  double SK6091::functionTest::sumG(Eigen::RowVector3d X, int n) {
	double f = 0.0, d = 0.0;
	for (size_t i = 1; i <=n; i++)
	{
		d += std::pow(CF(i, X), 2.0);
	}
	f = d;
	return f;
}
inline Eigen::Matrix<double, 2, 1>SK6091::functionTest::fDerp(Eigen::Matrix<double, 1, 2> X, Eigen::Matrix<double, 1, 5> y_i, Eigen::Matrix<double, 1, 5> t_i) {
	Eigen::Matrix<double, 2, 5> gradF;
	for (size_t i = 0; i < 1; i++)
	{
		for (size_t j = 0; j < 5; j++) {
			gradF(i, j) = std::exp(X(0, 1) * t_i(i, j));
		}
	}
	for (size_t i = 1; i < 2; i++)
	{
		for (size_t j = 0; j < 5; j++) {
			gradF(i, j) = X(0,0)*t_i(i-1,j)*std::exp(X(0, 1) * t_i(i-1, j));
		}
	}
	Eigen::Matrix<double, 5, 1> F;
	for (size_t i = 0; i < 5; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			F(i, j) = X(0,0)*std::exp(X(0, 1) * t_i(0, i)) -y_i(0,i);
		}
	}
	auto sum = gradF * F;
	return sum;
}
inline Eigen::Matrix<double, 2, 2> SK6091::functionTest::hes(Eigen::Matrix<double, 1, 2> X, Eigen::Matrix<double, 1, 5> y_i, Eigen::Matrix<double, 1, 5> t_i) {
	Eigen::Matrix<double, 2, 5> gradF;
	for (size_t i = 0; i < 1; i++)
	{
		for (size_t j = 0; j < 5; j++) {
			gradF(i, j) = std::exp(X(0, 1) * t_i(i, j));
		}
	}
	for (size_t i = 1; i < 2; i++)
	{
		for (size_t j = 0; j < 5; j++) {
			gradF(i, j) = X(0, 0) * t_i(i - 1, j) * std::exp(X(0, 1) * t_i(i - 1, j));
		}
	}
	auto sum = gradF * gradF.transpose();
	return sum;
}
inline Eigen::Matrix<double, 5, 1>SK6091::functionTest::H_x(Eigen::Matrix<double, 1, 2> X, Eigen::Matrix<double, 1, 5>y_i, Eigen::Matrix<double, 1, 5> t_i) {
	Eigen::Matrix<double, 5, 1> temp;
	for (size_t i = 0; i < 5; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			temp(i, j) = X(0, 0)  * std::exp(X(0, 1) * t_i(0, i))  - y_i(0,i);
		}
	}
	return temp;
}
//test for 11 data
inline Eigen::Matrix<double, 4, 1>SK6091::functionTest::fDerp(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
	Eigen::Matrix<double, 4, 11> gradF;
	for (size_t i = 0; i < 1; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = std::exp(X(0, 1) * t_i(i, j));
		}
	}
	for (size_t i = 1; i < 2; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = t_i(0,j)* X(0, 0)* std::exp(X(0, 1) * t_i(0, j));
		}
	}

	for (size_t i = 2; i < 3; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = std::exp(X(0, 3) * t_i(0, j));
		}
	}
	for (size_t i = 3; i < 4; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = t_i(0, j) * X(0, 2) * std::exp(X(0, 3) * t_i(0, j));
		}
	}


	Eigen::Matrix<double, 11, 1> F;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			F(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i))  + X(0, 2) * std::exp(X(0, 3) * t_i(0, i))
				- y_i(0, i);
		}
	}
	auto sum = gradF * F;
	return sum;
}
inline  Eigen::Matrix<double, 4, 4> SK6091::functionTest::hes(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
	Eigen::Matrix<double, 4, 11> gradF;
	for (size_t i = 0; i < 1; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = std::exp(X(0, 1) * t_i(i, j));
		}
	}
	for (size_t i = 1; i < 2; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = t_i(0, j) * X(0, 0) * std::exp(X(0, 1) * t_i(0, j));
		}
	}

	for (size_t i = 2; i < 3; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = std::exp(X(0, 3) * t_i(0, j));
		}
	}
	for (size_t i = 3; i < 4; i++)
	{
		for (size_t j = 0; j < 11; j++) {
			gradF(i, j) = t_i(0, j) * X(0, 2) * std::exp(X(0, 3) * t_i(0, j));
		}
	}
	auto sum = gradF * gradF.transpose();
	return sum;

}


#endif	