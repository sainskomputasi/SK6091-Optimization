#pragma once
#ifndef SPIRAL_HEADER
#define SPIRAL_HEADER
#include <cmath>
#include <math.h>
#include <eigen/dense>
#include "spiralImp.h"
namespace SK6091 {
	class Spiral {
	public:
		inline static double rosenBrock(Eigen::RowVector3d, int);
		inline static double t_rosenBrock(Eigen::RowVector2d,  int );
		inline static double constrainedFun(Eigen::RowVector2d, int);

		inline static double f_rosenBrock(Eigen::RowVector4d, int);


		static Eigen::RowVector3d center(Eigen::Matrix<double, 50, 3>, int);
		static Eigen::RowVector3d move(Eigen::RowVector3d,Eigen::Matrix<double,3,3>, Eigen::Matrix<double, 3, 3>, Eigen::RowVector3d);
		static Eigen::RowVector3d solve(int);

		static Eigen::RowVector2d t_center(Eigen::Matrix<double, 50, 2>, int);
		static Eigen::RowVector2d t_move(Eigen::RowVector2d, Eigen::Matrix<double, 2, 2>, Eigen::Matrix<double, 2, 2>, Eigen::RowVector2d);
		static Eigen::RowVector2d t_solve(unsigned int,int);

		inline static double LeastSquareModel(Eigen::RowVector2d x, Eigen::Matrix<double, 1, 5>, Eigen::Matrix<double, 1, 5>);
		static Eigen::RowVector2d t_center(Eigen::Matrix<double, 50, 2>, int, Eigen::Matrix<double, 1, 5>, Eigen::Matrix<double, 1, 5>);
		static Eigen::RowVector2d t_solve(unsigned int, Eigen::Matrix<double, 1, 5>, Eigen::Matrix<double, 1, 5>);
		
		//in r4 
		static Eigen::RowVector4d f_center(Eigen::Matrix<double, 50, 4>,   int);
		static Eigen::RowVector4d f_move(Eigen::RowVector4d, Eigen::Matrix<double, 4, 4>, Eigen::Matrix<double, 4, 4>, Eigen::RowVector4d);
		static Eigen::RowVector4d f_solve(int);
		//nonlinear least squares 
		inline static double fLeastSquares(Eigen::RowVector4d x, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i);
		static Eigen::RowVector4d l_center(Eigen::Matrix<double, 250, 4>, int, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i);
		static Eigen::RowVector4d l_solve(unsigned int, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i);


	};
}
#endif