#ifndef SPECIAL_FUNCTION
#define SPECIAL_FUNCTION
#include <cmath>
#include <eigen/dense>
namespace SK6091 {
	class functionTest
	{
	public:
		double Rosenbrock(Eigen::RowVector2d );
		double Griewank(Eigen::RowVector2d);
		double Ackley(Eigen::RowVector2d);
		double Booth_s(Eigen::RowVector2d);
		double Michalewicz(Eigen::RowVector2d);
		double Spring_sys(Eigen::RowVector2d);
		double Rastrigin(Eigen::RowVector2d);


	private:
	
	};


}
#include "specialFunctionImp.hpp"
#endif