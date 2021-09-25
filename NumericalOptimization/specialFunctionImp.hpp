#ifndef SPECIAL_FUNCTION_IMP
#define SPECIAL_FUNCTION_IMP
#include <cmath>
#include <math.h>
#include <eigen/dense>
#include "specialFunction.hpp"
double SK6091::functionTest::Rosenbrock(Eigen::RowVector2d temp) {
	return (1+ (1/4000) *std::pow(temp[0],2.0) + (1.0 / 4000.0) *std::pow(temp[1], 2.0) - std::cos(temp[0]) * std::cos((1.0/2.0) * temp[1]*std::sqrt(2.0)));
}
double SK6091::functionTest::Griewank(Eigen::RowVector2d temp){
	return (100*(std::pow(temp[1]-std::pow(temp[0],2.0),2.0)) + std::pow(1-temp[0],2.0)); 
}
double SK6091::functionTest::Ackley(Eigen::RowVector2d temp) {
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
double SK6091::functionTest::Booth_s(Eigen::RowVector2d temp) {
	auto s = 0.0;
	s = std::pow((temp[0] + 2 * temp[1] - 7.0), 2.0) + std::pow((2 * temp[0] + temp[1] - 5.0),2.0);
	return s;
}
double SK6091::functionTest::Michalewicz(Eigen::RowVector2d temp) {
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
double SK6091::functionTest::Spring_sys(Eigen::RowVector2d temp) {
	return (100 * std::pow(std::sqrt(std::pow(temp[0], 2.0) + std::pow((temp[1] + 1.0), 2.0)) - 1.0, 2.0) + 90 * std::pow(std::sqrt(std::pow(temp[0], 2.0) + std::pow((temp[1] - 1.0), 2.0)) - 1.0, 2.0) - (20 * temp[0] + 40 * temp[1]));
}
double SK6091::functionTest::Rastrigin(Eigen::RowVector2d temp) {
	auto begin = 0, end = 2;
	auto s = 0.0;
	while (begin!=end)
	{
		s = s + (std::pow(temp[begin], 2.0) - 10 * std::cos(2 * M_PI * temp[begin]) );
		++begin;
	}
	return (10*2+s);
}
#endif