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
		inline static Eigen::Matrix<double, 1, 4>Levenberg(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double,1, 11>);
		inline static Eigen::Matrix<double, 4, 1> levenbergJacobian(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::Matrix<double, 4, 4>levenbergHessian(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static double fx(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		//declaration for dongleg 
		inline static Eigen::Matrix<double, 1, 4>dogLeg(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::Matrix<double, 4, 1> dogLegJacobian(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::Matrix<double, 4, 4>dogLegHessian(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static double dogLegfx(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::Matrix<double, 11, 1>dogLegFx(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::Matrix<double, 4, 11> dogLegtempJacob(Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 1, 11>, Eigen::Matrix<double, 1, 11>);
		inline static Eigen::RowVector3d homotopy(Eigen::RowVector3d);

		inline static double Ft(int, Eigen::RowVector2d);
		inline static double Pt(int, int, Eigen::RowVector2d);
		inline static Eigen::RowVector2d homotopyt(Eigen::RowVector2d);

		inline static Eigen::RowVector2d polakRibere(Eigen::RowVector2d);
		//implementation of levenberg method.. 
		inline static Eigen::RowVector2d levenvergMarquardt(Eigen::RowVector2d);
		inline static Eigen::Matrix<double, 1, 4> gaussNewton(Eigen::Matrix<double, 1, 4> X, Eigen::Matrix<double, 1, 11> y_i, Eigen::Matrix<double, 1, 11>t_i);

	private:
		std::string notFound;
	};


}
#include "specialFunctionImp.hpp"
#endif