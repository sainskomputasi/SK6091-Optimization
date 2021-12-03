#ifndef MULTI_DIMEN
#define MULTI_DIMEN
#include <iostream>
#include<vector>
#include <cmath>
#include <Eigen/dense>
namespace SK6091{
    class MultiD{
        private:

        public:
            Eigen::RowVector2d stepestDes(Eigen::RowVector2d);
            Eigen::RowVector2d newton(Eigen::RowVector2d);
            Eigen::RowVector2d newton(Eigen::RowVector2d,int);
            Eigen::RowVector2d quasiNewton(Eigen::RowVector2d);
            inline double backtrackingLineSearch(Eigen::RowVector2d, Eigen::RowVector2d);
    };
}
#include "multiDimenImp.hpp"
#endif