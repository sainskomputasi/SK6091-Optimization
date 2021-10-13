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
				term1 = Griewank(temp);
				temp[i] = vecName[i] - del;
				term2 = Griewank(temp);
				term3 = Griewank(vecName);
				temporary(i, j) = (term1 - 2.0 * term3 + term2) / (std::pow(del, 2.0));
			}
			else {
				temp = vecName;
				temp[i] = vecName[i] + del;
				temp[j] = vecName[j] + del;
				term1 = Griewank(temp);

				temp = vecName;
				temp[i] = vecName[i] + del;
				temp[j] = vecName[j] - del;
				term2 = Griewank(temp);

				temp = vecName;
				temp[i] = vecName[i] - del;
				temp[j] = vecName[j] + del;
				term3 = Griewank(temp);


				temp = vecName;
				temp[i] = vecName[i] - del;
				temp[j] = vecName[j] - del;
				term4 = Griewank(temp);
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


	Eigen::Matrix<double, 11, 1> F;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			F(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i))
				- y_i(0, i);
		}
	}
	auto sum = gradF * F;
	return sum;
}

//test for 11 data
inline Eigen::Matrix<double, 4, 1>SK6091::functionTest::levenbergJacobian(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
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

//levenberg implementation 
inline Eigen::Matrix<double, 4, 1> levenbergJacobian(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11>y_i, Eigen::Matrix<double, 1, 11> t_i) {
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


	Eigen::Matrix<double, 11, 1> F;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			F(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i))
				- y_i(0, i);
		}
	}
	auto sum = gradF * F;
	return sum;
}

inline  Eigen::Matrix<double, 4, 4>SK6091::functionTest::levenbergHessian(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11>y_i, Eigen::Matrix<double, 1, 11>t_i) {
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
inline double SK6091::functionTest::fx(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
	Eigen::Matrix<double, 11, 1> temp;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			temp(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) - y_i(0, i);
		}
	}
	double sigma = 0.0;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++)
		{
			sigma += std::pow(temp(i, j), 2.0);

		}
	}
	return ((1 / 2.0) * sigma);
}

inline Eigen::Matrix<double, 1, 4>SK6091::functionTest::Levenberg(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
	//levenberg implementation 
	unsigned int k = 0;
	double v = 2.0;
	Eigen::Matrix<double, 1, 4> xNew;
	Eigen::Matrix<double, 4, 4> A;
	A = SK6091::functionTest::levenbergHessian(X, y_i, t_i);
	Eigen::Matrix<double, 4, 1>g;
	g = SK6091::functionTest::levenbergJacobian(X, y_i, t_i);
	const double tolerance = 1.e-8;
	bool found = (g.lpNorm<Eigen::Infinity>()<=tolerance);
	const double  tau = 1.e-9;
	Eigen::Matrix<double, 1, 4> maxTemp;
	maxTemp << A(0, 0), A(1, 1), A(2, 2), A(3, 3);
	double mu = tau * A.maxCoeff();
	double gainRatio = 0.0,gainRationTerm=0.0;
	Eigen::Matrix<double, 4, 1> hlm;
	Eigen::Matrix<double, 4, 4> identity;
	identity = identity.Identity();
	double maxTerm = 0.0;
	while ((!found && k<100))
	{
		++k;
		std::cout << "total iteration \t:" << k << std::endl;
		hlm = (A + mu * A.Identity()).inverse() * -g;
		if ((hlm.norm()< (tolerance*(X.norm()+tolerance)) ))
		{
			found = true;
			break;
		}
		else{
			xNew = X + hlm.transpose();
			gainRationTerm = (1 / 2.0) * hlm.transpose() * (mu * hlm - g);
			gainRatio = (SK6091::functionTest::fx(X, y_i, t_i) - SK6091::functionTest::fx(xNew, y_i, t_i)) / (gainRationTerm);
			if (gainRatio>0)
			{
				X = xNew;
				A = SK6091::functionTest::levenbergHessian(X, y_i, t_i);
				g = SK6091::functionTest::levenbergJacobian(X, y_i, t_i);
				found = (g.lpNorm<Eigen::Infinity>() < tolerance);
				maxTerm = 1 - std::pow(2.0 * gainRatio - 1.0, 3.0);
				mu = mu * std::max(1.0 / 3.0, maxTerm);
				v = 2.0;
			}
			else {
				mu = mu * v;
				v = 2.0 * v;
			}
		}
		std::cout << "test for x ....";
		std::cout << X;
	}
	return X;
}

//dog Leg implementation 
inline Eigen::Matrix<double, 4, 1>SK6091::functionTest::dogLegJacobian(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11>y_i, Eigen::Matrix<double, 1, 11> t_i) {
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


	Eigen::Matrix<double, 11, 1> F;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			F(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i))
				- y_i(0, i);
		}
	}
	auto sum = gradF * F;
	return sum;
}

inline  Eigen::Matrix<double, 4, 4>SK6091::functionTest::dogLegHessian(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11>y_i, Eigen::Matrix<double, 1, 11>t_i) {
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
inline double SK6091::functionTest::dogLegfx(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
	Eigen::Matrix<double, 11, 1> temp;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			temp(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) - y_i(0, i);
		}
	}
	double sigma = 0.0;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++)
		{
			sigma += std::pow(temp(i, j), 2.0);

		}
	}
	return (1 / 2.0 * sigma);
}

inline Eigen::Matrix<double, 1, 4>SK6091::functionTest::dogLeg(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11> t_i) {
	//dogleg implementation
	unsigned int k = 0;
	double delta = 1.0;
	Eigen::Matrix<double, 1, 4> xNew;
	Eigen::Matrix<double, 4, 1>g,hsd,hdl,hgn,a,b;
	g = SK6091::functionTest::dogLegJacobian(X, y_i, t_i);
	double alpha = 0.0,c=0.0;
	double tolerance = 1.e-6;
	double gainRatio = 0.0, gainRationTerm = 0.0;

	bool found = (SK6091::functionTest::dogLegFx(X, y_i, t_i).lpNorm<Eigen::Infinity>()<tolerance) || (g.lpNorm<Eigen::Infinity>()<tolerance) ;
	double beta = 0.0;
	while (!found && (k<200))
	{
		std::cout << "dogleg iteration  \t: " << k << std::endl;
		++k;
		alpha = std::pow(g.norm(), 2.0) / 
			((SK6091::functionTest::dogLegtempJacob(X, y_i, t_i).transpose() * g).norm());
		hsd = -alpha * g;
		hgn = SK6091::functionTest::dogLegtempJacob(X, y_i, t_i).inverse() * 
				-SK6091::functionTest::dogLegFx(X, y_i, t_i);
		a = alpha * hsd;
		b = hgn;
		c = a.transpose() * (b - a);
		if (c<0)
		{
			beta = (-c + std::sqrt(std::pow(c,2.0) + std::pow((b-a).norm(),2.0) *(std::pow(delta,2.0) - std::pow(a.norm(),2.0)) ) ) / 
				(std::pow((b-a).norm(),2.0));
		}
		else {
			beta = (std::pow(delta,2.0)  - std::pow(a.norm(),2.0))/
				(c + std::sqrt(std::pow(c, 2.0) + std::pow((b - a).norm(), 2.0) * (std::pow(delta, 2.0) - std::pow(a.norm(), 2.0))));
		}
		if (hgn.norm()<-delta)
		{
			hdl = hgn;
		}
		else if ((alpha * hsd.norm()) >= delta) {
			hdl = (delta / hsd.norm()) * hsd;
		}
		else {
			hdl = alpha * hsd + beta * (hgn - alpha * hsd);
		}
		if (hdl.norm() <(tolerance*(X.norm()+ tolerance)))
		{
			found = true;
		}
		else {
			xNew = X + hdl.transpose();
			//gain ration .... 
			if (hdl == hgn)
			{
				gainRationTerm = SK6091::functionTest::dogLegfx(X, y_i, t_i);
			}
			else if (hdl == (-delta * g / g.norm())) {
				gainRationTerm = delta * (2.0 * alpha * g.norm() - delta)/ (2.0*alpha);
			}
			else {
				gainRationTerm = (1 / 2.0 * alpha) * std::pow((1 - beta), 2.0) * std::pow(g.norm(), 2.0) + 
					beta * (2.0 - beta) * SK6091::functionTest::dogLegfx(X, y_i, t_i);
			}
			gainRatio = (SK6091::functionTest::dogLegfx(X, y_i, t_i) - SK6091::functionTest::dogLegfx(xNew, y_i, t_i))/
				gainRationTerm;
			if (gainRatio  >0 )
			{
				X = xNew;
				g = SK6091::functionTest::dogLegJacobian(X, y_i, t_i);
				found = (SK6091::functionTest::dogLegFx(X, y_i, t_i).lpNorm<Eigen::Infinity>()  < tolerance ) 
					|| (g.lpNorm<Eigen::Infinity>()<tolerance);
			}
			if (gainRatio>0.75)
			{
				delta = std::max(delta, 3.0 * hdl.norm());
			}
			else if (gainRatio < 0.25) {
				delta = delta / 2.0;
				found = (delta <= (tolerance * (X.norm() + tolerance)));
			}
		}

	}
	return X;
}
inline  Eigen::Matrix<double, 11, 1>SK6091::functionTest::dogLegFx(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11>y_i, Eigen::Matrix<double, 1, 11>t_i) {

	Eigen::Matrix<double, 11, 1> F;
	for (size_t i = 0; i < 11; i++)
	{
		for (size_t j = 0; j < 1; j++) {
			F(i, j) = X(0, 0) * std::exp(X(0, 1) * t_i(0, i)) + X(0, 2) * std::exp(X(0, 3) * t_i(0, i))
				- y_i(0, i);
		}
	}
	return F;
}
inline Eigen::Matrix<double, 4, 11>SK6091::functionTest::dogLegtempJacob(Eigen::Matrix<double, 1, 4>X, Eigen::Matrix<double, 1, 11>y_i, Eigen::Matrix<double, 1, 11> t_i) {
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

	return gradF;
}

inline Eigen::RowVector3d SK6091::functionTest::homotopy(Eigen::RowVector3d X1) {
	const int N = 4;
	const double h = 1.0 / 4.0;  //step length for runga kutta.. 
	unsigned int iter = 1, bound = 4,kk=1;
	double RMM[4];
	Eigen::RowVector3d x, b,ki;
	Eigen::Matrix<double, 3, 3> A,Ainv;
	Eigen::Matrix<double, 4, 10>xk1;
	RMM[0] = 0.5;
	RMM[1] = 0.5;
	RMM[2] = 1.0;
	RMM[3] = 0.0;
	for (size_t i = 1; i < N; i++)
	{
		x[i - 1] = X1[i - 1]; // move the initial value to temporary object x..
	}
	for (size_t i = 1; i <N	; i++)
	{
		b[i - 1] = -h*SK6091::functionTest::F(i, x);  //step length *f(x)
	}

	while (iter<=bound) // make 4*4 system of linear eq.. 
	{
		for (size_t i = 1; i <N; i++)
		{
			x[i - 1] = X1[i - 1]; //update value for x ... contain row vector 1*3
		}
		kk = 1;
		while (kk<=4)
		{
			for (size_t i = 1; i <N; i++)
			{
				for (size_t j = 1; j <N; j++) {
					A(i - 1, j - 1) = SK6091::functionTest::P(i, j, x);//jacobian matrix... 
				}
			}
			//solve ax  = b here ... 
			Ainv = A.inverse();//note for singular... assume that A is singular matrix 
			ki = (Ainv * b.transpose()).transpose(); //find each value for k_i

			for (size_t i = 1; i <N; i++)
			{
				xk1(kk-1,i-1) = ki[i - 1];
				x[i - 1] = X1[i-1]+RMM[kk-1]*xk1(kk-1,i-1); 
			}
			++kk;//increase kk and move to next rows
		}
		for (size_t i = 1; i <N; i++)
		{
			X1[i-1]  = X1[i - 1] + (xk1(0, i - 1) + 2 * xk1(1, i - 1) + 2 * xk1(2, i - 1) + xk1(3,i-1)) / 6.0; //apply runga kutta orde 4
		}
		
			std::cout << "x" << std::endl;
			std::cout<< X1 << std::endl;
		
		//update solution....
		++iter;
	}
	return X1;
}
//homotopy case for n = 2
inline double SK6091::functionTest::Ft(int I, Eigen::RowVector2d X) {
	double f = 0.0;

	switch (I) {
	case 1:
		f = std::pow(X[0], 2.0) + std::pow(X[1], 2.0) - 1;
		break;
	case 2:
		f = std::sin(X[0]) - X[1]; //original function
		break;

	}
	return f;

}
inline  double SK6091::functionTest::Pt(int I, int J, Eigen::RowVector2d X) {
	double p = 0.0;

	switch (I) {
		case 1:
			switch (J) {
			case 1:
				p = 2*X[0];
				break;
			case 2:
				p = 2*X[1];
				break;
			}
			break;
		case 2:
			switch (J) {
				case 1:
					p = std::cos(X[0]);
					break;
				case 2:
					p = -1;
					break;
			}
			break;
	}
	return p;
}

inline Eigen::RowVector2d SK6091::functionTest::homotopyt(Eigen::RowVector2d X1) {
	const int N = 3;
	const double h = 1.0 / 4.0;  //step length for runga kutta.. 
	unsigned int iter = 1, bound = 4, kk = 1;
	double RMM[4];
	Eigen::RowVector2d x, b, ki;
	Eigen::Matrix<double, 2, 2> A, Ainv;
	Eigen::Matrix<double, 4, 10>xk1;
	RMM[0] = 0.5;
	RMM[1] = 0.5;
	RMM[2] = 1.0;
	RMM[3] = 0.0;
	for (size_t i = 1; i < N; i++)
	{
		x[i - 1] = X1[i - 1]; // move the initial value to temporary object x..
	}
	for (size_t i = 1; i < N; i++)
	{
		b[i - 1] = -h * SK6091::functionTest::Ft(i, x);  //step length *f(x)
	}

	while (iter <= bound) // make 4*4 system of linear eq.. 
	{
		for (size_t i = 1; i < N; i++)
		{
			x[i - 1] = X1[i - 1]; //update value for x ... contain row vector 1*3
		}
		kk = 1;
		while (kk <= 4)//design for rungga kutta method..
		{
			for (size_t i = 1; i < N; i++)
			{
				for (size_t j = 1; j < N; j++) {
					A(i - 1, j - 1) = SK6091::functionTest::Pt(i, j, x);//jacobian matrix... 
				}
			}
			//solve ax  = b here ... 
			Ainv = A.inverse();//note for singular... assume that A is singular matrix 
			ki = (Ainv * b.transpose()).transpose(); //find each value for k_i

			for (size_t i = 1; i < N; i++)
			{
				xk1(kk - 1, i - 1) = ki[i - 1];
				x[i - 1] = X1[i - 1] + RMM[kk - 1] * xk1(kk - 1, i - 1);
			}
			++kk;//increase kk and move to next rows
		}
		for (size_t i = 1; i < N; i++)
		{
			X1[i - 1] = X1[i - 1] + (xk1(0, i - 1) + 2 * xk1(1, i - 1) + 2 * xk1(2, i - 1) + xk1(3, i - 1)) / 6.0; //apply runga kutta orde 4
		}

		std::cout << "x" << std::endl;
		std::cout << X1 << std::endl;

		//update solution....
		++iter;
	}
	return X1;
}
//conjugate gradient method
//polak robiere implementation
//case Griewank function -- for d =2;
inline Eigen::RowVector2d SK6091::functionTest::polakRibere(Eigen::RowVector2d X) {
	const double tolerance = 1.e-7;
	auto begin = 1, end = 100;//max iteration
	double falPrev = SK6091::functionTest::Griewank(X);
	//APPLY ALGORITHM 5.4 - case for ploak ribere ..
	Eigen::RowVector2d deriv,derivPrev,searchPrev,search;
	double alpha, falpha;
	//benchmark purpose
	std::ofstream polakAnalis;
	polakAnalis.open("polakTable.csv", std::ios::app);
	while (begin!=end)
	{

		if (begin==1)
		{
			derivPrev = SK6091::functionTest::grad(X);//aprocimation for gradient
			searchPrev = -derivPrev;
			alpha = SK6091::functionTest::goldFunc(X, searchPrev)[0];//find alpha using line search method.. 
			falpha = SK6091::functionTest::goldFunc(X, searchPrev)[1];
			if (derivPrev.norm()<tolerance)
			{
				return X;// for k =1.. step if satisfy condition (see termination condition nocedal book )
				break;//first guess.. stop
			}
			X += alpha * searchPrev; //update x for next iteration.
			falPrev = SK6091::functionTest::Griewank(X);

		}
		else {
			deriv = SK6091::functionTest::grad(X);//aprimation for gradient when x is updated
			//evaluate f'(k+1) (5.41a -Nocedal book)
			//change bkfr to polak ribere method bkPR (5.4.1a - Nocedal book, page. 52)
			search = -deriv + (std::pow(deriv.norm(),2.0) / 
				(std::pow(derivPrev.norm(), 2.0))) * searchPrev;//update pk+1, searchPrev == pk
			//point 5.41b nocedal book (algorithm 5.4 FR)
			//termination condition. 
			if (deriv.norm() < tolerance)
			{
				return X;
				break;
			}
			//else, update the current x to new X 
			derivPrev = deriv;
			searchPrev = search;//note compound statement x= x+alpha*search
			X += alpha * search;  // update  for next X... and return x*
		}
		polakAnalis <<begin<<";"<< X[0] << ";" << X[1] << std::endl;
		++begin;//update k (5.4.1c - nocedal book)
	}
	polakAnalis.close();
	return	X;
}
//definition of levenberg method.. 
inline Eigen::RowVector2d SK6091::functionTest::levenvergMarquardt(Eigen::RowVector2d X) {
	double mu = 1.0	;
	const double tolerance = 1.e-18;
	unsigned int begin = 1, end = 50;
	Eigen::Matrix<double, 2, 2> secDerph,identity;
	Eigen::RowVector2d search,derph;
	double r = 0.0;
	double fx = 0.0, fPrev = 0.0;
	bool posDefinite = false;
	while (begin!=end)
	{

		fPrev = SK6091::functionTest::Griewank(X);
		derph = SK6091::functionTest::grad(X);
		secDerph = SK6091::functionTest::hessian(X);//aproximation for hessian, gradient.
		
		while (!(secDerph + mu * secDerph.Identity()).determinant() > 0)
		{
			std::cout << "increase mu" << std::endl;
			mu = 2.0 * mu ;
		}
		search = ((secDerph + mu * secDerph.Identity()).inverse() * -derph.transpose()).transpose();
		r = (fPrev - SK6091::functionTest::Griewank(X+search))  /
			( fPrev - (fPrev+search*derph.transpose() +1/2.0*search*secDerph*search.transpose()));
		if (r>1.e-3)
		{
			X = X+search;
			mu = mu * std::max(1.0 / 3.0, 1.0 - std::pow(2 * r - 1.0, 3.0)); //f decrease
		}
		else {
			mu = 2.0 * mu; //increase mu but keep x
		}
		std::cout << "total iteration \t: " << begin << "\t mu\t: " << mu <<"\t x[0]\t:" <<X[0]<<"\t x[1]\t: " <<X[1] << std::endl;
		//termination condition.. 
		if (derph.lpNorm<Eigen::Infinity>() <tolerance)
		{
			return X;
			break;
		}
		++begin;
	}
	return X;
}


#endif	