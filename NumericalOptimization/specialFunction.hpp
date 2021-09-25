#ifndef SPECIAL_FUNCTION
#define SPECIAL_FUNCTION
#include <cmath>
#include <eigen/dense>
namespace SK6091 {
	class functionTest
	{
	public:
		inline double Rosenbrock(Eigen::RowVector2d );
		inline double Griewank(Eigen::RowVector2d);
		inline double Ackley(Eigen::RowVector2d);
		inline double Booth_s(Eigen::RowVector2d);
		inline double Michalewicz(Eigen::RowVector2d);
		inline double Spring_sys(Eigen::RowVector2d);
		inline double Rastrigin(Eigen::RowVector2d);


	private:
	
	};


}
#include "specialFunctionImp.hpp"
#endif