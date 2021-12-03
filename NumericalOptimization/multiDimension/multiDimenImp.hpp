#ifndef MULTI_DIMEN_IMP
#define MULTI_DIMEN_IMP
#include <iostream>
#include<vector>
#include <cmath>
#include "multiDimen.hpp"
#include "../specialFunctionImp.hpp"
#include <fstream>
Eigen::RowVector2d SK6091::MultiD::stepestDes(Eigen::RowVector2d x){
    double eps=1e-6,eps1=eps;
    double falpPrev = SK6091::functionTest::Griewank(x);
    int begin=1;
    const int MAX=3000;
    Eigen::RowVector2d deriv,search;
    deriv << 0, 0;
    search << 0, 0;
    auto tempAlpha = deriv, tempFalph = deriv;
    double alpha1 = 0.0, falpha1 = 0.0;
    std::ofstream write; //file handling 
    write.open("DataAnalisisTri.csv", std::ios::app);
    write << "Iteration;x1;x2;f(x);Norm" << std::endl;

    while (begin!=MAX)
    {
        deriv = SK6091::functionTest::grad(x);
        search = -deriv;
        tempAlpha = SK6091::functionTest::goldFunc(x, search);
        alpha1 = tempAlpha[0];
        falpha1 = tempAlpha[1];
        write << begin << ";" << x[0] << ";" << x[1] << ";" << falpha1 << ";" << deriv.norm() << ";" << std::endl;
        if ((std::fabs(falpha1 - falpPrev) < eps) || (deriv.norm() < eps1)) {
            return x;
            break;
        }

        falpPrev = falpha1;
        x = x + alpha1 * search;
        ++begin;
    }
    

    return x ;


}
Eigen::RowVector2d SK6091::MultiD::newton(Eigen::RowVector2d x) {
    auto eps = 1e-7;
    double fPrev =0.0,f= 0.0;
    auto deriv = x;
    deriv << 0, 0;
    Eigen::Matrix2d hes,inv;
    hes << 0, 0, 0, 0;
    inv << 0, 0, 0, 0;
    int begin = 1,end=500;
    Eigen::Vector2d notCoverge;
    while (begin!=end)
    {
        fPrev = SK6091::functionTest::Griewank(x);
        deriv = SK6091::functionTest::grad(x);
        hes = SK6091::functionTest::hessian(x);
        inv = hes.inverse();
        x = (x.transpose() - (inv * deriv.transpose())).transpose();
        f = SK6091::functionTest::Griewank(x);
        std::cout << "i :" << begin << "\tx : " << x << "\tf: " << f << "\tnorm : " << deriv.norm() << std::endl;
        if ((std::fabs(f-fPrev) < eps) || (deriv.norm()<eps))
        {
            return x;
            break;
        }
        ++begin;
    }
    if (begin==(end-1))
    {
        std::cout << "Method fail to coverge ...." << std::endl;
    }
    return notCoverge;
}
Eigen::RowVector2d SK6091::MultiD::newton(Eigen::RowVector2d x, int maxIter) {
    auto eps = 1e-3;
    double fPrev = 0.0, f = 0.0;
    auto deriv = x;
    deriv << 0, 0;
    Eigen::Matrix2d hes, inv;
    hes << 0, 0, 0, 0;
    inv << 0, 0, 0, 0;
    int begin = 1;
    auto search = x;
    search << 0, 0;
    Eigen::Vector2d notCoverge;
    double alpha1=0.0, falpha1 = 0.0;
    std::ofstream write; //file handling 
    write.open("DataAnalisisTri.csv", std::ios::app);
    write << "Iteration;x1;x2;f(x);Norm" << std::endl;
    while (begin != maxIter)
    {
        fPrev = SK6091::functionTest::Griewank(x);
        deriv = SK6091::functionTest::grad(x);
        hes = SK6091::functionTest::hessian(x);
        inv = hes.inverse();
        search = -inv * deriv.transpose();
        std::cout << "backracking start .... " << std::endl;
        alpha1 = SK6091::MultiD::backtrackingLineSearch(x,search);
        //alpha1 = SK6091::functionTest::goldFunc(x, search)[0];

        std::cout << "backracking finish .... with alpha \t:  " <<alpha1<< std::endl;

        //x = (x.transpose() - (inv * deriv.transpose())).transpose();
        f = SK6091::functionTest::Griewank(x);
        std::cout << begin << ";" << x[0] << ";" << x[1] << ";" << falpha1 << ";" << deriv.norm() << ";" << std::endl;
        if ( (deriv.norm() < eps))
        {
            return x;
            std::cout << "Total iteration \t: " << begin << "result \t: " << x << std::endl;
            break;
        }
        fPrev = falpha1;
        x = x + alpha1 * search;
        falpha1 = SK6091::functionTest::Griewank(x);
        ++begin;
    }
    if (begin == (maxIter - 1))
    {
        std::cout << "Method fail to coverge ...." << std::endl;
    }
    return notCoverge;
}
Eigen::RowVector2d SK6091::MultiD::quasiNewton(Eigen::RowVector2d x) {
    auto tolerance = 1e-6;
    Eigen::Matrix2d A;
    A = A.Identity();
    auto fPrev = SK6091::functionTest::Griewank(x);
    auto alpha1 = 0.0, falpa = 0.0;
    auto derivPrev = x;
    auto search = x; 
    auto begin = 1, end = 100;
    auto deltax = x;
    auto deriv = x;
    auto deltag = x;
    auto term1 = A, term2 = A;
    std::ofstream write; //file handling 
    write.open("DataAnalisisTri.csv", std::ios::app);
    write << "Iteration;x1;x2;f(x);Norm" << std::endl;
    //write << "x;y;iteration" << std::endl;
    while (begin!=end)
    {
        if (begin == 1) {
            derivPrev = SK6091::functionTest::grad(x);
            search = -derivPrev;//steepestDes
            alpha1 = SK6091::functionTest::goldFunc(x, search)[0];
            falpa = SK6091::functionTest::goldFunc(x, search)[1];
            if (std::fabs(falpa - fPrev) < 1e-3) {
                std::cout << "phase 1\t: pass ..." << std::endl;
                return x;
                break;
            }
            fPrev = falpa;
            x = x + alpha1 * search;
        }
        else { //case for 2 .. flow here
            deltax = (alpha1 * search);
            if (begin > 2) {
                deltax = deltax.transpose();
                search = search.transpose();
            }
            //debug for 2nd iteration.. test...
            deriv = SK6091::functionTest::grad(x);
            deltag = deriv - derivPrev;//notice for transf ...
            term1 = (deltag.transpose() * deltag) / (deltag*deltax.transpose());
            term2 = (derivPrev.transpose() * derivPrev) / (derivPrev * search.transpose());
            A = A+term1 + term2;//debug here, consider matrix A
            search = -(A.inverse()) * deriv.transpose();//pass
            alpha1 = SK6091::functionTest::goldFunc(x, search.transpose())[0];
            falpa = SK6091::functionTest::goldFunc(x, search.transpose())[1];
            write << begin<<";"<<x[0] << ";" << x[1] << ";" << falpa << ";" << deriv.norm() << ";" << std::endl;
            //std::cout << "Falpha \t: " << falpa << "fPrev\t: " << fPrev << std::endl;
            if ((std::fabs(falpa - fPrev) < tolerance) || (deriv.norm() < tolerance)) {
                return x;
                break;
            }
            fPrev = falpa;
            derivPrev = deriv;
            x += alpha1 * search;
        }
            ++begin;
    }
    write.close();
    return x;
}
inline double SK6091::MultiD::backtrackingLineSearch(Eigen::RowVector2d x, Eigen::RowVector2d s) {
    double alpha = 1.0;
    double p = 0.1;
    double c =0.9;
    int begin = 0;
    while (SK6091::functionTest::Griewank(x+alpha*s) >= 
        (SK6091::functionTest::Griewank(x) + c*alpha*SK6091::functionTest::grad(x)*s.transpose()) )
    {
        ++begin;
        std::cout << "++begin\t: " << ++begin << std::endl;
        alpha=alpha*p;
    }
    return alpha;
}

#endif