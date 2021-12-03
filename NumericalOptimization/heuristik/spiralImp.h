#pragma once
#ifndef SPIRAL_IMPLEMENTATION
#define SPIRAL_IMPLEMENTATION
#include <cmath>
#include <math.h>
#include <eigen/dense>
#include "spiralHeader.h"
#include <iostream>
#include "../specialFunctionImp.hpp"
#include <random>
inline double SK6091::Spiral::rosenBrock(Eigen::RowVector3d x, int n) {
	double z = 0.0;
	double a=0.0, b = 0.0;
	for (size_t i = 0; i < n-1; i++)
	{
		a = 100.0 * (std::pow(( x[i+1]  - std::pow(x[i],2.0)),2.0));
		b = std::pow((1 - x[i]), 2.0);
		z += a + b;
	}
	double err = 0.0;
	err = std::pow((z - 0.0), 2.0);
	return err;
}
inline double SK6091::Spiral::t_rosenBrock(Eigen::RowVector2d x,  int n) {
	double z = 0.0;
	double a = 0.0, b = 0.0;
	for (size_t i = 0; i < n - 1; i++)
	{
		a = 100.0 * (std::pow((x[i + 1] - std::pow(x[i], 2.0)), 2.0));
		b = std::pow((1 - x[i]), 2.0);
		z += a + b;
	}  //griewank test....
	//double z=(1 + (1 / 4000) * std::pow(x[0], 2.0) + (1.0 / 4000.0) * std::pow(x[1], 2.0) - std::cos(x[0]) * std::cos((1.0 / 2.0) * x[1] * std::sqrt(2.0)));
	double err = 0.0;
	err = std::pow((z - 0.0), 2.0);
	return err;
}

inline  double SK6091::Spiral::constrainedFun(Eigen::RowVector2d x, int n) {
	/*double fx = std::pow(x[0], 2.0) + std::pow(x[1], 2.0);
	double constraintE = std::pow(x[0] - 1.0, 3.0) - std::pow(x[1], 2.0);
	double alpha = 1.e10;
	fx = fx + alpha * std::pow(constraintE, 2.0); */ //ex :1 x* (1,0)
	double beta = 1.e20;
	
	double fx = 1 / 2.0 * std::pow(x[0] - 2.0, 2.0)+1 / 2.0 * std::pow(x[1] - 1 / 2.0, 2.0);
	double c1 = -std::pow(x[0] + 1.0, -1.0) + x[1] + 1 / 4.0;
	double c2 = -x[0];
	double c3 = -x[1];  //example 1
	/*fx = 2.0 * (std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 1.0) - x[0];
	c1 = std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 1.0;*/ //example 2

	//fx = x[0] + x[1];
	//c1 = std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 2.0;
	fx = std::pow(x[0] - 3.0, 2.0) + std::pow(x[1] - 1.0, 2.0); //problem 3
	c1 = (std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 16.0);
	c2 = (4.0 - std::pow(x[0], 2.0) - std::pow(x[1], 2.0));
	c3 = 0;

	fx = std::pow(x[0] - 3.0, 2.0) + std::pow(x[1] - 1.0, 2.0); //problem 1
	c1 = 1.0-std::pow(x[0],2.0)-std::pow(x[1]-2.0,2.0);
	c2 = x[1]-4.0;
	c3 = -x[1] + (1.0 / 4.0) * std::pow(x[0], 2.0);

	//problem 5
	fx = std::pow(x[0], 2.0) + std::pow(x[1] - 1.0, 2.0);
	c1 = x[1] - std::pow(x[0], 2.0);
	double sumConstrained = std::pow(c1, 2.0);
	//problem 2
	fx = std::pow(x[0] - 10.0, 2.0) + std::pow(x[1] - 20.0, 2.0);
	c1 = -std::pow(x[0] - 5.0, 2.0) - std::pow(x[1] - 5.0, 2.0) + 100.0;
	c2 = std::pow(x[0] - 6.0, 2.0) + std::pow(x[1] - 5.0, 2.0) - 82.81;
	//double 
	sumConstrained = beta * std::pow(std::max(0.0, c1), 2.0) + beta * std::pow(std::max(0.0, c2), 2.0);	//beta * std::pow(std::max(0.0, c2), 2.0)+
		//; //inequality constrained (penalty) 
	//problem nocedal page 441prob 6
	fx = 2.0 * (std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 1.0) - x[0];
	c1 = std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 1.0;
	sumConstrained = beta*std::fabs(c1);
	//problem  7 igor griva 
	fx = std::exp(2.0 * x[0]) + std::exp(-4.0 * x[1]);
	c1 = std::pow(x[0], 2.0) + std::pow(x[1], 2.0) - 1.0;
	sumConstrained = beta * std::fabs(c1);
	fx = fx + sumConstrained;

	

	double err = 0.0;
	err = std::pow((fx - 0.0), 2.0);
	return err;
}

Eigen::RowVector3d SK6091::Spiral::center(Eigen::Matrix<double, 50, 3> x, int m) {
	int n = 3;
	double bestErr = -SK6091::functionTest::gSum(x.row(0));
	int idx = 0;
	double err = 0.0;
	for (int i = 0; i < m; i++)
	{
		err = -SK6091::functionTest::gSum(x.row(i));
		if (err<bestErr)
		{
			idx = i;
			bestErr = err;
		}
	}
	return x.row(idx);
}
Eigen::RowVector3d SK6091::Spiral::move(Eigen::RowVector3d x, Eigen::Matrix<double, 3, 3> r, Eigen::Matrix<double, 3, 3> rIdentity, Eigen::RowVector3d center) {
	Eigen::RowVector3d temp1,temp2;
	temp1 = (rIdentity * center.transpose()).transpose();
	temp2 = (r * x.transpose() - temp1.transpose()).transpose();
	return temp2;
}

Eigen::RowVector3d SK6091::Spiral::solve(int maxIteration) {
	double theta = M_PI/4.0, r = 0.94;
	int m = 50;
	int n = 3;
	Eigen::Matrix<double, 3, 3>R12, R13, R23;
	R12 << std::cos(theta), -std::sin(theta), 0,
		std::sin(theta), std::cos(theta), 0,
		0, 0, 1;
	R13 << std::cos(theta), 0, -std::sin(theta),
		0, 1, 0,
		std::sin(theta), 0, std::cos(theta);
	R23 << 1, 0, 0,
		0, std::cos(theta), -std::sin(theta),
		0, std::sin(theta), std::cos(theta);
	Eigen::Matrix<double, 3, 3> Rm;
	Rm = R23 * (R13 * R12);
	Rm = r * Rm;
	Eigen::Matrix<double, 3, 3> rMin;
	rMin = Rm - Rm.Identity();
	Eigen::Matrix<double, 50, 3> points;
	double lower_bound = -40.0;
	double upper_bound = 40.0;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
	double a_random_double = unif(re);
	//points = points.setRandom();
	for (int i = 0; i < points.rows(); i++)
	{
		for (int j = 0; j < points.cols(); j++)
		{
			a_random_double = unif(re);
			points(i, j) = a_random_double;
		}
	}
	Eigen::RowVector3d center;
	center = SK6091::Spiral::center(points, m);
	std::cout << "test best center \t: " << center << std::endl;
	Eigen::RowVector3d x;
	x << 0, 0, 0;
	for (size_t i = 0; i < maxIteration; i++)
	{
		for (size_t j = 0; j < m; j++) {
			x = points.row(j);
			points.row(j) = SK6091::Spiral::move(x, Rm, rMin, center);
		}
		center = SK6091::Spiral::center(points, m);
		std::cout << "test \t: " << SK6091::functionTest::gSum(center)<<std::endl;
	}
	return center;
}
//2 dimension 
inline Eigen::RowVector2d SK6091::Spiral::t_center(Eigen::Matrix<double, 50, 2>x, int m) {
	int n = 2;
	double bestErr = SK6091::Spiral::constrainedFun(x.row(0),n);
	int idx = 0;
	double err = 0.0;
	for (int i = 0; i < m; i++)
	{
		err = SK6091::Spiral::constrainedFun(x.row(i),2);
		if (err < bestErr)
		{
			idx = i;
			bestErr = err;
		}
	}
	return x.row(idx);
}
inline Eigen::RowVector2d SK6091::Spiral::t_move(Eigen::RowVector2d x, Eigen::Matrix<double, 2, 2>r, Eigen::Matrix<double, 2, 2>rIdentity, Eigen::RowVector2d center) {
	Eigen::RowVector2d temp1, temp2;
	temp1 = (rIdentity * center.transpose()).transpose();
	temp2 = (r * x.transpose() - temp1.transpose()).transpose();
	return temp2;
}
inline Eigen::RowVector2d SK6091::Spiral::t_solve(unsigned int maxIteration,int q) {
	double theta = M_PI / 4.0, r = 0.98;
	int m = 50;
	int n = 2;
	Eigen::Matrix<double, 2, 2>R12;
	R12 << std::cos(theta), -std::sin(theta),
		std::sin(theta), std::cos(theta);
	R12 = r * R12;
	Eigen::Matrix<double, 2, 2> rMin;
	rMin = R12 - R12.Identity();
	Eigen::Matrix<double, 50, 2> points;
	double lower_bound = -0.8;
	double upper_bound =0.8;
	srand(time(NULL));
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
	unif.reset();
	double a_random_double = unif(re);
	//points = points.setRandom();
	for (int i = 0; i < points.rows(); i++)
	{
		for (int j = 0; j < points.cols(); j++)
		{
			a_random_double = unif(re);
			points(i, j) = a_random_double;
		}
	}
	Eigen::RowVector2d center;
	center = SK6091::Spiral::t_center(points, m);
	std::cout << "test best center \t: " << center << std::endl;
	Eigen::RowVector2d x;
	x << 0, 0;
	for (size_t i = 0; i < maxIteration; i++)
	{
		for (size_t j = 0; j < m; j++) {
			x = points.row(j);
			points.row(j) = SK6091::Spiral::t_move(x, R12, rMin, center);
		}
		center = SK6091::Spiral::t_center(points, m);
		std::cout << "test \t: " << std::sqrt(SK6091::Spiral::constrainedFun(center,2)) << std::endl;
		//std::cout << "test \t: " << std::sqrt(SK6091::functionTest::t_gSum(center)) << std::endl;

	}
	return center;
}
//test least squares...
inline  double SK6091::Spiral::LeastSquareModel(Eigen::RowVector2d x, Eigen::Matrix<double, 1, 5> t_i,
	Eigen::Matrix<double, 1, 5>y_i) {
	double sum = 0.0;
	
	for (size_t i = 0; i < 5; i++)
	{
		//sum += std::pow(x[0] * std::exp(x[1] * t_i(0, i))  -y_i(0, i),2.0);
		//another model 
		sum += std::pow(x[0]*std::cos(x[1]*t_i(0,i)) + x[1]*std::sin(x[0]*t_i(0,i)) - y_i(0, i), 2.0);

	}
	sum = 1.0 / (1.0 + sum);
	double err = 0.0;
	err = std::pow((sum - 0.0), 2.0);
	return err;

}

inline Eigen::RowVector2d SK6091::Spiral::t_center(Eigen::Matrix<double, 50, 2>x, int m, Eigen::Matrix<double, 1, 5> t_i, Eigen::Matrix<double, 1, 5> y_i) {
	int n = 2;
	double bestErr = -SK6091::Spiral::LeastSquareModel(x.row(0),t_i,y_i);
	int idx = 0;
	double err = 0.0;
	for (int i = 0; i < m; i++)
	{
		err = -SK6091::Spiral::LeastSquareModel(x.row(i), t_i,y_i);
		if (err < bestErr)
		{
			idx = i;
			bestErr = err;
		}
	}
	return x.row(idx);
}

inline Eigen::RowVector2d SK6091::Spiral::t_solve(unsigned int maxIteration, Eigen::Matrix<double, 1, 5> t_i, Eigen::Matrix<double, 1, 5> y_i) {
	double theta = M_PI / 4.0, r = 0.94;
	int m = 50;
	int n = 2;
	Eigen::Matrix<double, 2, 2>R12;
	R12 << std::cos(theta), -std::sin(theta),
		std::sin(theta), std::cos(theta);
	R12 = r * R12;
	Eigen::Matrix<double, 2, 2> rMin;
	rMin = R12 - R12.Identity();
	Eigen::Matrix<double, 50, 2> points;
	double lower_bound = 0.0;
	double upper_bound = 110.0;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
	double a_random_double = unif(re);
	//points = points.setRandom();
	for (int i = 0; i < points.rows(); i++)
	{
		for (int j = 0; j < points.cols(); j++)
		{
			a_random_double = unif(re);
			points(i, j) = a_random_double;
		}
	}
	Eigen::RowVector2d center;
	center = SK6091::Spiral::t_center(points, m,t_i,y_i);
	std::cout << "test best center \t: " << center << std::endl;
	Eigen::RowVector2d x;
	x << 0, 0;
	for (size_t i = 0; i < maxIteration; i++)
	{
		for (size_t j = 0; j < m; j++) {
			x = points.row(j);
			points.row(j) = SK6091::Spiral::t_move(x, R12, rMin, center);
		}
		center = SK6091::Spiral::t_center(points, m, t_i, y_i);
		//std::cout << "test \t: " << std::sqrt(SK6091::Spiral::t_rosenBrock(center,2)) << std::endl;
		//std::cout << "test \t: " << std::sqrt(SK6091::functionTest::t_gSum(center)) << std::endl;
		//std::cout << "value \t: " << std::sqrt(SK6091::Spiral::LeastSquareModel(center, t_i, y_i)) << std::endl;;

	}
	return center;
}
//test for r-4 
inline double SK6091::Spiral::f_rosenBrock(Eigen::RowVector4d x, int n) {
	double z = 0.0;
	double a = 0.0, b = 0.0;
	for (size_t i = 0; i < (n - 1); i++)
	{
		a = 100.0 * (std::pow((x[i + 1] - std::pow(x[i], 2.0)), 2.0));
		b = std::pow((1 - x[i]), 2.0);
		z += a + b;
	}  //griewank test....
	//double z=(1 + (1 / 4000) * std::pow(x[0], 2.0) + (1.0 / 4000.0) * std::pow(x[1], 2.0) - std::cos(x[0]) * std::cos((1.0 / 2.0) * x[1] * std::sqrt(2.0)));
	double err = 0.0;
	err = std::pow((z - 0.0), 2.0);
	return err;
	//grewank test ;
	/*unsigned int begin = 0, end = n - 1;
	double s = 0.0;
	while (begin != end)
	{
		s += 1 + (1 / 4000) * std::pow(x[begin], 2.0) + (1.0 / 4000.0) * std::pow(x[begin + 1], 2.0) - std::cos(x[begin]) * std::cos((1.0 / 2.0) * x[begin + 1] * std::sqrt(2.0));
		++begin;
	}
	double err = 0.0;
	err = std::pow((s - 0.0), 2.0);
	return err;
	*/
}
Eigen::RowVector4d SK6091::Spiral::f_center(Eigen::Matrix<double, 50, 4> x,  int m) {
	int n = 4;
	double bestErr = SK6091::Spiral::f_rosenBrock(x.row(0), n);
	int idx = 0;
	double err = 0.0;
	for (int i = 0; i < m; i++)
	{
		err = SK6091::Spiral::f_rosenBrock(x.row(i), n);
		if (err < bestErr)
		{
			idx = i;
			bestErr = err;
		}
	}
	return x.row(idx);
}

Eigen::RowVector4d SK6091::Spiral::f_move(Eigen::RowVector4d x, Eigen::Matrix<double, 4, 4> r, Eigen::Matrix<double, 4, 4>rIdentity, Eigen::RowVector4d center) {
	
	Eigen::RowVector4d temp1, temp2;
	temp1 = (rIdentity * center.transpose()).transpose();
	temp2 = (r * x.transpose() - temp1.transpose()).transpose();
	return temp2;
 }
Eigen::RowVector4d SK6091::Spiral::f_solve(int maxIteration) {
	double theta = M_PI / 4.0, r = 0.98;
	int m = 50;
	int n = 4;

	Eigen::Matrix<double, 4, 4>R12,R13,R14,R23,R24,R34,R43;

	R12 <<	std::cos(theta), -std::sin(theta), 0, 0,
			std::sin(theta), std::cos(theta), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;
	R13 << std::cos(theta), 0, -std::sin(theta), 0,
		0, 1, 0, 0,
		std::sin(theta), 0, std::cos(theta), 0,
		0, 0, 0, 1;
	R14 << std::cos(theta), 0, 0, -std::sin(theta),
		0, 1, 0, 0,
		0, 0, 1, 0,
	std::sin(theta), 0, 0, std::cos(theta);

	R23 << 1, 0, 0, 0,
		0, std::cos(theta), -std::sin(theta), 0,
		0, std::sin(theta), std::cos(theta), 0,
		0, 0, 0, 1;

	R24 << 1, 0, 0, 0,
		0, std::cos(theta), 0, -std::sin(theta),
		0, 0,1, 0,
		0, std::sin(theta), 0, std::cos(theta);

	R34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0,	0, std::cos(theta), -std::sin(theta),
		0, 0, std::sin(theta), std::cos(theta);

	R13 = r * (R34 * (R24 * (R23 * (R14 * (R13 * R12)))));//store in r13

	Eigen::Matrix<double, 4, 4> rMin;
	rMin = R13 - R13.Identity();
	Eigen::Matrix<double, 50, 4> points;
	double lower_bound = -2.0;
	double upper_bound = 2.0;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
	double a_random_double = unif(re);
	//points = points.setRandom();
	for (int i = 0; i < points.rows(); i++)
	{
		//a_random_double = unif(re);

		for (int j = 0; j < points.cols(); j++)
		{
			a_random_double = unif(re);
			points(i, j) = a_random_double;
		}
		std::cout<<std::endl;
	}
	Eigen::RowVector4d center;
	center = SK6091::Spiral::f_center(points, m);
	std::cout << "test best center \t: " << center << std::endl;
	Eigen::RowVector4d x;
	x << 0, 0,0,0;
	for (size_t i = 0; i < maxIteration; i++)
	{
		for (size_t j = 0; j < m; j++) {
			x = points.row(j);
			points.row(j) = SK6091::Spiral::f_move(x, R13, rMin, center);
		}
		std::cout << std::sqrt(SK6091::Spiral::f_rosenBrock(center, 4)) << "  " << std::endl;

		center = SK6091::Spiral::f_center(points, m);
		std::cout<< std::endl;
	}
	std::cout << "result point \t: " << center << std::endl;
	return center;
}

//implement least squares r-4
inline  double SK6091::Spiral::fLeastSquares(Eigen::RowVector4d x, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i) {
	double sum = 0.0;
	for (size_t i = 0; i < 11; i++)
	{
		sum += 1/2.0 *std::pow(x[0] * std::exp(x[1] * t_i(0, i)) + x[2]*std::exp(x[3]*t_i(0,i)) -y_i(0,i), 2.0);
	}
	//sum = 1.0 / (1.0 + sum);
	//sum = 1 / 2.0 * sum;
	double err = 0.0;
	err = std::pow((sum - 0.0), 2.0);
	return err;
}

Eigen::RowVector4d SK6091::Spiral::l_center(Eigen::Matrix<double, 250, 4>x, int m, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i) {
	int n = 4;
	double bestErr = SK6091::Spiral::fLeastSquares(x.row(0), y_i, t_i);
	int idx = 0;
	double err = 0.0;
	for (int i = 0; i < m; i++)
	{
		err = SK6091::Spiral::fLeastSquares(x.row(i), y_i, t_i);
		if (err < bestErr)
		{
			idx = i;
			bestErr = err;
		}
	}
	return x.row(idx);
}
Eigen::RowVector4d SK6091::Spiral::l_solve(unsigned int maxIteration, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i) {
	double theta = M_PI / 4.0, r = 0.98;
	int m = 250;
	int n = 4;

	Eigen::Matrix<double, 4, 4>R12, R13, R14, R23, R24, R34, R43;

	R12 << std::cos(theta), -std::sin(theta), 0, 0,
		std::sin(theta), std::cos(theta), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;
	R13 << std::cos(theta), 0, -std::sin(theta), 0,
		0, 1, 0, 0,
		std::sin(theta), 0, std::cos(theta), 0,
		0, 0, 0, 1;
	R14 << std::cos(theta), 0, 0, -std::sin(theta),
		0, 1, 0, 0,
		0, 0, 1, 0,
		std::sin(theta), 0, 0, std::cos(theta);

	R23 << 1, 0, 0, 0,
		0, std::cos(theta), -std::sin(theta), 0,
		0, std::sin(theta), std::cos(theta), 0,
		0, 0, 0, 1;

	R24 << 1, 0, 0, 0,
		0, std::cos(theta), 0, -std::sin(theta),
		0, 0, 1, 0,
		0, std::sin(theta), 0, std::cos(theta);

	R34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, std::cos(theta), -std::sin(theta),
		0, 0, std::sin(theta), std::cos(theta);

	R13 = r * (R34 * (R24 * (R23 * (R14 * (R13 * R12)))));//store in r13

	Eigen::Matrix<double, 4, 4> rMin;
	rMin = R13 - R13.Identity();
	Eigen::Matrix<double, 250, 4> points;
	double lower_bound = -1.0;
	double upper_bound = 1.0;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;
	double a_random_double = unif(re);
	//points = points.setRandom();
	for (int i = 0; i < points.rows(); i++)
	{
		//a_random_double = unif(re);

		for (int j = 0; j < points.cols(); j++)
		{
			a_random_double = unif(re);
			points(i, j) = a_random_double;
		}
		std::cout << std::endl;
	}
	Eigen::RowVector4d center;
	center = SK6091::Spiral::l_center(points, m,y_i,t_i);
	std::cout << "test best center \t: " << center << std::endl;
	Eigen::RowVector4d x;
	x << 0, 0, 0, 0;
	for (size_t i = 0; i < maxIteration; i++)
	{
		for (size_t j = 0; j < m; j++) {
			x = points.row(j);
			points.row(j) = SK6091::Spiral::f_move(x, R13, rMin, center);
		}
		std::cout << std::sqrt(SK6091::Spiral::fLeastSquares(center,y_i,t_i)) << "  center \t: " <<center << std::endl;

		center = SK6091::Spiral::l_center(points, m,y_i,t_i);
		
		std::cout << std::endl;
	}
	std::cout << "result point \t: " << center << std::endl;
	return center;
}

#endif