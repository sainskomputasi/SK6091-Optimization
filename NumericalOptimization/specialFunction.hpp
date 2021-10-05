#ifndef SPECIAL_FUNCTION
#define SPECIAL_FUNCTION
#include <cmath>
#include <eigen/dense>
#include <fstream>
namespace SK6091 {
	class functionTest
	{
	public:
		inline  static double Rosenbrock(Eigen::RowVector2d );
		inline static double Griewank(Eigen::RowVector2d);
		inline static double Ackley(Eigen::RowVector2d);
		inline static double Booth_s(Eigen::RowVector2d);
		inline static double Michalewicz(Eigen::RowVector2d);
		inline static double Spring_sys(Eigen::RowVector2d);
		inline static double Rastrigin(Eigen::RowVector2d);
		inline static Eigen::RowVector2d grad(Eigen::RowVector2d vecName);
		static Eigen::Matrix2d  hessian(Eigen::RowVector2d vecName);
		static Eigen::RowVector2d goldFunc(Eigen::RowVector2d, Eigen::RowVector2d);
		inline static double F(int, Eigen::RowVector3d);
		inline static double P(int, int, Eigen::RowVector3d);
		inline static double CF(int, Eigen::RowVector3d);
		inline static double  sumDerph(int, Eigen::RowVector3d);
		inline static double sumG(Eigen::RowVector3d, int);
		//test case for least squares.... 
		inline static Eigen::Matrix<double,2,1> fDerp(Eigen::Matrix<double, 1, 2>, Eigen::Matrix<double,1,5>, Eigen::Matrix<double, 1, 5>);
		inline static Eigen::Matrix<double, 2, 2> hes(Eigen::Matrix<double, 1, 2>, Eigen::Matrix<double, 1, 5>, Eigen::Matrix<double, 1, 5>);
		inline static Eigen::Matrix<double, 5, 1>H_x(Eigen::Matrix<double, 1, 2>,Eigen::Matrix<double,1,5>, Eigen::Matrix<double, 1, 5>);
		inline static Eigen::Matrix<double, 4, 1> fDerp(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::Matrix<double, 4, 4> hes(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);

	private:
		std::string notFound;
	};


}
#include "specialFunctionImp.hpp"
#endif