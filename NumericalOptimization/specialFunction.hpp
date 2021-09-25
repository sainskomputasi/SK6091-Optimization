#ifndef SPECIAL_FUNCTION
#define SPECIAL_FUNCTION
#include <cmath>
#include <eigen/dense>
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

	private:
		std::string notFound;
	};


}
#include "specialFunctionImp.hpp"
#endif