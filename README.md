# SK6091-Optimization
  **SK6091-Optim** is a lightweight C++ library of numerical optimization methods and part of the paper (**Finding Solution Of Nonlinear Constrained Optimization Using Spiral Optimization Algorithm**), it is written by Michael sipayung. The goals are to complete an independent research project and to build an Optimization library for C++ which can provide solutions for a single variable and multi-variable optimization. The users should be able to select the algorithm to be used and also access any class abstraction. Available Algorithms: golden section search, dichotomous Search, Fibonacci Search, Steepest-Descent Method, Newton Method, Gauss-Newton Method, quasi newton method, Particle Swarm Optimization, Simulated annealing, Spiral Dynamic Optimization, Differential Evolution, Evolutionary algorithms. I also provide a benchmark for each algorithm at the end of each section that consists of the optimization algorithm. The following graph is intended to be used as a sample to compile. 
  
  ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/ezgif.com-gif-maker%20(6).gif)
  

## 1.The optimization Problem
  - **Introduction**

      Optimization problems of sorts arise in all quantitative disciplines from computer science and engineering to operations research and economics. In the simplest case, an optimization problem consists of maximizing or minimizing a real function by systematically choosing input values from within an allowed set and computing the value of the function. More generally, optimization includes finding "best available" values of some objective function given a defined domain (or input), including a variety of different types of objective functions and different types of domains. The following graph is illustrated local and global minima/ maxima.
        ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/dsBuffer.bmp.png)

  - **The Basic Optimization Problem & C++** 

      When using this library/ module, you start by writting an objective function that may be single variable or many variable, this objective function always give a return value and normally it's return double data type, double data type has a good preccision rather than float and sometimes it's size equals in the memory. The folowing code is a simple way to return a value from an objective function **(rastrigin function)** using c++. 
      ```c++
      double oneDrastrigin(double x) {
	    return (10+std::pow(x,2.0)-10*std::cos(2*M_PI*x));
      }
      ```
      ```c++
      double twoDrastrigin(double x,double y) {
	    return (20+ (std::pow(x,2.0)-10*std::cos(2*M_PI*x)+
            std::pow(y,2.0)-10*std::cos(2*M_PI*y)));
      }
      ```
      Rather than defining a new function for each type, we can define a function template. A function template is formula from which we can generate type-specific version of that functions. The template version of Rastarign function looks like

      ```c++
      template <typename T>
      inline double oneD_Rastrign(const T& x) {
	    return (10 + std::pow(x, 2.0) - 10 * std::cos(2 * M_PI * x));
      }
      ```
      ```c++
      template <typename T>
      inline double twoD_Rastrign(const T& x,const T& y) {
	    return (20+ (std::pow(x,2.0)-10*std::cos(2*M_PI*x)+
                std::pow(y,2.0)-10*std::cos(2*M_PI*y)));
      }
      ```
      **sk6091::optimization** library/ module also uses the preprocessor to define header guards. Header guards rely on preprocessor variables. Preprocessor variables have one of two possible states: define or not defined. sk6091::optimization use these facilitiesto guard against multiple inclussion as follows:
      ```c++        
    #ifndef ONE_DIMENSIONAL
    #define ONE_DIMENSIONAL
    #include <cmath>
    #define MAGN(a,b) ((b)>=0.0?std::fabs(a):-std::fabs(a))
    #define GOLD 0.3819660
    #define ZEPS 1.0e-10
    #define SHIFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
    namespace SK6091 {
	    class OneDimension {
	    public:
		    double brent(double(*)(double), double, double, double,double ,double *, int);

	    private:
	    };
    }
    #include "oneDimensional.hpp"
    #endif
       ```
  - **How to compile**

    Prerequisites
    - Windows 7 or newer (Linux)
    - Git
    - Visual Studio or Gcc

    First download sk6091::optimization itself, it can be installed anywhere. We recommend somewhere like `C:\optim\testOptimization` or `D:\optim\testOptimization` 
    
    ```cmd
    > git clone https://github.com/sainskomputasi/SK6091-Optimization.git
    ```
    Now let's make a simple CMake project with a main file.
    ```cmake
    # CMakeLists.txt
    cmake_minimum_required(VERSION 3.0)
    project(test)

    add_executable(optimization NumericalOptimization.cpp)
    ```
    ```c++
    //NumericalOptimization.cpp
    #include "oneDimen/oneDimensional.hpp"
    double objectiveFunc(double x) {
	    return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0);
    }
    int main()
    {
	    auto xmi = -1.5;
	    SK6091::OneDimension testBrent;
	    std::cout << "minimum occurs at x \t: " << 
            testBrent.brent(objectiveFunc, -4.0, -2.0, 0.01, 0.0001, &xmi, 100);
    }
     ```
    Then, we build our project in the normal CMake way:
    ```cmd
    > cmake ..
    ```



## 2.One Dimensional Optimization
  - **Introduction**
    
    This section will describe techniques to find the minimum or maximum of a function of a single variable, f(x).The one-dimensional (1-D) optimization problem refers to an objective function with one variable. In practice, optimization problems with many variables are complex, and rarely does one find a problem with a single variable. However, 1-D optimization algorithms form the basic building blocks for 
    multivariable algorithms. As these algorithms form a subproblem of multivariable optimization problems, numerous methods (or algorithms) have been reported in the literature, each with some unique advantage over the others. As an example, a single-variable objective function could be 
    > **f(x) = 2*sin(x)-((x^2)/10.0)** 
    ; Using Gnu-Plot to plot this function will look like 
    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/plot.jpg)
    
 
    The following benchmark results have been generated using some numerical experiments.
    
    N   | Method Name   | Maximum Point | Minimum Point |Tolerance |Iteration|Time       
    ----|---------------|---------------|---------------|----------|---------|-------------
    1|  Newton|1.42755| *|1.00E-04|
    2|  Secant|1.42755| *|1.00E-04|
    3|  Golden Section|1.42755|-1.74633|1.00E-04|(35,35)|0.1831 s
    4|  Parabolic Interpolation|1.42755| *|1.00E-04 
    5|  Brent Method| *|-1.74633|1.00E-04|(8)|0.1978 s
    6|  Newton|1.42755| *| 1E-11|
    7|  Secant|1.42755| *| 1E-11|
    8|  Golden Section|1.42755|-1.74633|1E-11|(63,64)|0.0257 s
    9|  Parabolic Interpolation|1.42755|*|1E-11| 
    10| Brent Method| *|-1.74633|1E-11|(10)|0.0474 s

  - **Newton Method**
    - Newton’s method is an open method similar to Newton-Raphson because it does not require initial guesses that bracket the optimum [a,b]. In the newton method in this library, we need not a first derivative and second derivative to calculate, because we already provide an approximation to them. See The finite-difference derivative approximations for calculating derivatives in the calculating derivative section.The following code is intended to demonstrate how we call the  Newton method through this module. *(Note : namespace for this module is **SK6091** and folowed by it's class)*.   
        ```c++
        #include "oneDimen/oneDimensional.hpp"
        double objectiveFunc(double x) {
	        return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0);
        }
        int main()
        {
	        SK6091::OneDimension testNewton;
	        auto guess = 0.5;
	        std::cout << "maximum occurs at x \t: " << 
                testNewton.newton(objectiveFunc, guess, 0.00001, 100);
	        return 0;
        }
        ```
        After compile the program : 
        > The maximum occurs at x : 1.42755

        **The Following code is Newton abstraction**
        ```c++
        #ifndef ONE_DIMENSIONAL
        #define ONE_DIMENSIONAL
        #include <cmath>
        #define R ((pow(5,0.5)-1.0)/2.0)
        namespace SK6091 {
	        class OneDimension {
	        public:
		        double newton(double(double),double &,double,int);
	        private:
	        };
        }
        #include "oneDimensional.hpp"
        #endif
        ```
    **The Newton function is consist of 4 parameters:**
    * *Objective Function*
    * *initial guess (reference type)*
    * *tolerance for error*
    * *max iteration*

  - **Brent's Method**
    - Brent's method combines the slow, dependable golden-section search with the faster, but possibly unreliable,parabolic interpolation.The algorithm tries to use the potentially fast-converging secant method or inverse quadratic interpolation if possible, but it falls back to the more robust bisection method if necessary. Matlab has a  professional minimization function that implements the advantages of brent's method called fminbnd. The following code is intended to demonstrate how we call the brent's  method through this module. *(Note : namespace for this module is **SK6091** and folowed by it's class)*.  
        ```c++
        #include "oneDimen/oneDimensional.hpp"
        double objectiveFunc(double x) {
	        return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0);
        }
        int main()
        {
	        auto xmi = -1.5;
	        SK6091::OneDimension testBrent;
	        std::cout << "minimum occurs at x \t: " << 
                testBrent.brent(objectiveFunc, -4.0, -2.0, 0.01, 0.0001, &xmi, 100);
        }
        ```
        After compile the program : 
        > The minimum occurs at x : -1.74633

        **The Following code is brent's method abstraction**
        ```c++
        #ifndef ONE_DIMENSIONAL
        #define ONE_DIMENSIONAL
        #include <cmath>
        #define MAGN(a,b) ((b)>=0.0?std::fabs(a):-std::fabs(a))
        #define GOLD 0.3819660
        #define ZEPS 1.0e-10
        #define SHIFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
        namespace SK6091 {
	        class OneDimension {
	        public:
		        double brent(double(*)(double), double, double, double,double ,double *, int);

	        private:
	        };
        }
        #include "oneDimensional.hpp"
        #endif
        ```
    **The brent function is consist of 6 parameters:**
    * *Objective Function (pointer to function type)*
    * *lower bound*
    * *auxiliary bound (f(lower)<f(aux)<f(upper))*
    * *upper bound*
    * *tolerance for error*
    * *initial guess (pointer type)*
    * *max iteration*
  - **Golden-Section Search**
    - The golden-section search is a simple, general-purpose, single-variable search technique. It is similar in spirit to the bisection approach for locating root. The following code is intended to demonstrate how we call the golden section method through this module. *(Note : namespace for this module is **SK6091** and folowed by it's class)*. 
        ```c++
        #include "oneDimen/oneDimensional.hpp"
        double objectiveFunc(double x) {
	        return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0);
        }
        int main(){
            SK6091::OneDimension testGolden;
	        std::cout<<"The maximum occurs at x \t: "<< 
                    testGolden.goldenSec(objectiveFunc, 0.0, 4.0, 100, 0.00001);
        }
        ```
        After compile the program : 
        > The maximum occurs at x : 1.42755

        **The Following code is golden section abstraction**
        ```c++
        #ifndef ONE_DIMENSIONAL
        #define ONE_DIMENSIONAL
        #include <cmath>
        #define R ((pow(5,0.5)-1.0)/2.0)
        namespace SK6091 {
	        class OneDimension {
	        public:
		        double goldenSec(double(double), double, double, int, double);
	        private:
	        };
        }
        #include "oneDimensional.hpp"
        #endif
        ```
    **The golden section function is consist of 5 parameters:**
    * *Objective Function*
    * *lower bound*
    * *upper bound*
    * *max iteration*
    * *tolerance for error*
  - **Secant Method**
    - Secant Method is a similiar to bisection method that require initial guesses that bracket the optimum [a,b]. The difference between secant method and the bisection method is in the secant method, both the magnitude and the sign of the derivative 
      are used to locate the zero of f′(x).The following code is intended to demonstrate how we can call the the secant method through this module/ library *(Note : namespace for this module is **SK6091** and folowed by it's class)* 
        ```c++
        #include "oneDimen/oneDimensional.hpp"
        double objectiveFunc(double x) {
	        return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0);
        }
        int main()
        {
	        SK6091::OneDimension testSecant;
	        std::cout << "maximum occurs at x \t: " << 
                testSecant.secant(objectiveFunc,0.0,4.0,0.00001, 100);
	        return 0;
        }
        ```
        After compile the program : 
        > The maximum occurs at x : 1.42755

        **The Following code is secant method abstraction**
        ```c++
        #ifndef ONE_DIMENSIONAL
        #define ONE_DIMENSIONAL
        #include <cmath>
        #define R ((pow(5,0.5)-1.0)/2.0)
        namespace SK6091 {
	        class OneDimension {
	        public:
		        double secant(double(double), double,double, double, int);
	        private:
	        };
        }
        #include "oneDimensional.hpp"
        #endif
        ```
    **The Secant function is consist of 5 parameters:**
    * *Objective Function*
    * *x0 as lower bound (double type)*
    * *x1 as higher bound (double type)*
    * *tolerance for error*
    * *max iteration*
  - **Parabolic Interpolation**
    - Parabolic interpolation takes advantage of the fact that a second-order polynomial often provides a good approximation to the shape of f(x) near an optimum. Parabolic interpolation will be converging rapidly if we care to select the best intervals. The following code is intended to demonstrate how we can call the parabolic interpolation method through this module *(Note : namespace for this module is **SK6091** and folowed by it's class)* 
        ```c++
        #include "oneDimen/oneDimensional.hpp"
        double objectiveFunc(double x) {
	        return (2 * std::sin(x) - std::pow(x, 2.0) / 10.0);
        }
        int main()
        {
	        SK6091::OneDimension testParabolic;
	        double x0 = 0.0, x1 = 1.0, x2 = 4.0;
	        std::cout << "maximum occurs at x \t: " << 
                        testParabolic.parabolic(objectiveFunc, x0, x1, x2, 100, 0.00001);
	        return 0;
        }
        ```
        After compile the program : 
        > The maximum occurs at x : 1.42755

        **The Following code is parabolic interpolation abstraction**
        ```c++
        #ifndef ONE_DIMENSIONAL
        #define ONE_DIMENSIONAL
        #include <cmath>
        #define R ((pow(5,0.5)-1.0)/2.0)
        namespace SK6091 {
	        class OneDimension {
	        public:
		        double parabolic(double(double), double&, double&, double&, int, double);
	        private:
	        };
        }
        #include "oneDimensional.hpp"
        #endif
        ```
    **The Parabolic interpolation function is consist of 6 parameters:**
    * *Objective Function*
    * *x0 as initial guesses (reference type)*
    * *x1 as initial guesses (reference type)*
    * *x2 as initial guesses (reference type)* 
    * *max iteration*
    * *tolerance for error*
## 3.Multi-Dimensional Gradient Method
  - **Introduction**

    As for one-dimensional optimization, there are two general classes of multidimensional methods,namely search method and gradient methods. The difference between gradient and direct search methods in multi-dimensional optimization is similar to the difference between these approaches in one-dimensional optimization. Direct search methods are useful when the derivative of the optimization function is not available to effectively guide the search for the optimum. While direct search methods explore the parameter space in a systematic manner, they are not computationally very efficient. On the other hand, gradient methods use information from the derivatives of the optimization function to more effectively guide the search and find optimum solutions much quicker. In this module/ library we also provided a saveral function that usually used in optimization. The Header specialFunction.hpp is The abstraction of the special function which is widely used in optimization and we choose to use inline for optimize the our code. 
    **The Following code is specialFunction.hpp implementation**

    ```c++
    #include "specialFunction.hpp"
    inline double SK6091::functionTest::Griewank(Eigen::RowVector2d temp){
	    return (1 + (1 / 4000) * std::pow(temp[0], 2.0) + (1.0 / 4000.0) * 
            std::pow(temp[1], 2.0) - std::cos(temp[0]) * std::cos((1.0 / 2.0) * 
                temp[1] * std::sqrt(2.0)));
    }
    inline  double SK6091::functionTest::Rosenbrock(Eigen::RowVector2d temp) {
	    return (100 * (std::pow(temp[1] - std::pow(temp[0], 2.0), 2.0)) + 
            std::pow(1 - temp[0], 2.0));
    }

    #endif
    ```
    The following benchmark results have been generated using some numerical experiments.

    
    Function   | Initial Point | Minimum Point |N | Initial Point | Minimum Point |N | Initial Point | Minimum Point |N |       
    ---------------|---------------|---------------|----------|---------|-------------|-------------|-------------|-------------|-------------
    Rosenbrock|1.9,1.9|0.999,0.998|10|0.5,0.5|0.999,0.998|27|1.1,1.2|1.000,1.000|21|
    Griewank|0.1,0.1|0.000,0.001|4|0.5,0.5|-0.00,0.00|4|1.1,1.2|0.0,-0.001|5|
    Ackley|0.1,-0.1|-0.00;0.000|4|0.3,0.3|-0.00;0.000|4|-0.2|-0.00;0.000|6|
    Booths|0.5,0.5|1.000,3.000|3|2.3,2.3|1.000,2.999|3|1.4, 2.0|0.999,3.000|4|
    Michalewicz|2.1,1.5|2.202,1.570|4|2.2,1.5|2.202,1.570|4|1.9,1.5|7.108,0.124|2|
    Rastrigin|0.1,-0.1|-0.99,0.99|2|0.0,-0.0|-1.3e-05|4|-0.2|-0.99,-0.99|2|
    Spring_sys|0.1,-0.1|0.504;0.121|8|0.0,-0.0|0.504;0.121|6|-0.2|0.504,0.121|8|

    The following graph give us knowledge how robust our module to estimate the optimum from special function which is widely used in optimization

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/benchmark%20special%20function.jpg)

  - **Steepest-Descent Method**
    - The steepest descent method ensures a reduction in the function value at every iteration. If the starting point is far away from the minimum, the gradient will be higher and the function reduction will be maximized in each iteration. The following code is intended to demonstrate how we can call the steepest descent method  to minimize Griewank function through this module *(Note : namespace for this module is **SK6091** and folowed by it's class)* 
        ```c++
        #include "multiDimension/multiDimenImp.hpp"
        int main()
        {
	        Eigen::RowVector2d initialPoint;
	        initialPoint << 0.5, 0.3;
	        SK6091::MultiD test;
	        std::cout<< "The Minimum Point occurs at x,y \t:\n"<<
                test.stepestDes(initialPoint)<<std::endl;
	        return 0;
        }
        ```
        After compile the program : 
        > The Minimum Point occurs at x,y :  1.02883 1.05855

        **The Following code is Steepest Descent abstraction**
        ```c++
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
            };
        }
        #include "multiDimenImp.hpp"
        #endif
        ```
        
        **The Following code is Griewank implementation**
        ```c++
        #ifndef SPECIAL_FUNCTION_IMP
        #define SPECIAL_FUNCTION_IMP
        #include <cmath>
        #include <math.h>
        #include <eigen/dense>
        #include "specialFunction.hpp"
        inline double SK6091::functionTest::Griewank(Eigen::RowVector2d temp){
	        return (1 + (1 / 4000) * std::pow(temp[0], 2.0) + (1.0 / 4000.0) * 
                std::pow(temp[1], 2.0) - std::cos(temp[0]) * std::cos((1.0 / 2.0) * 
                    temp[1] * std::sqrt(2.0)));
        }
        #endif
        ```     
        
  - **Newton Method**
      - Newton Method, the search direction in this method is based on the first and second derivative information. In this method, if the initial starting point is far away from the optimum, the search direction may not always be descent. Often a restart is required with a different starting point to avoid this difficulty. The following code is intended to demonstrate how we can call the Newton method  to minimize Griewank function through this module *(Note : namespace for this module is **SK6091** and folowed by it's class)* 
        ```c++
        #include "multiDimension/multiDimenImp.hpp"
        int main()
        {
	        Eigen::RowVector2d initialPoint;
	        initialPoint << 0.5, 0.3;
	        SK6091::MultiD test;
	        std::cout<< "The Minimum Point occurs at x,y \t:\n"<<
                test.newton(initialPoint)<<std::endl;
	        return 0;
        }
        ```
        After compile the program : 
        > The Minimum Point occurs at x,y :  0.9998 0.9996

        **The Following code is Newton abstraction**
        ```c++
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
                    Eigen::RowVector2d newton(Eigen::RowVector2d);
            };
        }
        #include "multiDimenImp.hpp"
        #endif
        ```
        
        **The Following code is Griewank implementation**
        ```c++
        #ifndef SPECIAL_FUNCTION_IMP
        #define SPECIAL_FUNCTION_IMP
        #include <cmath>
        #include <math.h>
        #include <eigen/dense>
        #include "specialFunction.hpp"
        inline double SK6091::functionTest::Griewank(Eigen::RowVector2d temp){
	        return (1 + (1 / 4000) * std::pow(temp[0], 2.0) + (1.0 / 4000.0) * 
                std::pow(temp[1], 2.0) - std::cos(temp[0]) * std::cos((1.0 / 2.0) * 
                    temp[1] * std::sqrt(2.0)));
        }
        #endif
        ```

  - **Quasi Newton Method**
     - Quasi-Newton Methods (QNMs) are generally a class of optimization methods that are used in Non-Linear Programming when full Newton’s Methods are either too time consuming or difficult to use. More specifically, these methods are used to find the global minimum of a function f(x) that is twice-differentiable. There are distinct advantages to using Quasi-Newton Methods over the full Newton's Method for expansive and complex non-linear problems. The following code is intended to demonstrate how we can call the Quasi Newton method  to minimize Griewank function through this module *(Note : namespace for this module is **SK6091** and folowed by it's class)*  
        ```c++
        #include "multiDimension/multiDimenImp.hpp"
        int main()
        {
	        Eigen::RowVector2d initialPoint;
	        initialPoint << 0.5, 0.3;
	        SK6091::MultiD test;
	        std::cout<< "The Minimum Point occurs at x,y \t:\n"<<
                test.quasiNewton(initialPoint)<<std::endl;
	        return 0;
        }
        ```
        After compile the program : 
        > The Minimum Point occurs at x,y :  1.00007 1.00018

        **The Following code is Quasi Newton abstraction**
        ```c++
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
                    Eigen::RowVector2d quasiNewton(Eigen::RowVector2d);
            };
        }
        #include "multiDimenImp.hpp"
        #endif
        ```
        
        **The Following code is Griewank implementation**
        ```c++
        #ifndef SPECIAL_FUNCTION_IMP
        #define SPECIAL_FUNCTION_IMP
        #include <cmath>
        #include <math.h>
        #include <eigen/dense>
        #include "specialFunction.hpp"
        inline double SK6091::functionTest::Griewank(Eigen::RowVector2d temp){
	        return (1 + (1 / 4000) * std::pow(temp[0], 2.0) + (1.0 / 4000.0) * 
                std::pow(temp[1], 2.0) - std::cos(temp[0]) * std::cos((1.0 / 2.0) * 
                    temp[1] * std::sqrt(2.0)));
        }
        #endif
        ```
## 4.Optimization - Nonlinear Least Squares 
  - **Introduction** 

    Non-linear least squares is the form of least squares analysis used to fit a set of m observations with a model that is non-linear in n unknown parameters (m ≥ n). It is used in some forms of nonlinear regression. The basis of the method is to approximate the model by a linear one and to refine the parameters by successive iterations. For example, we had collected data consisting of the size of a population of antelope of various times. Here t corresponds to the time at which the population y was counted. Suppose we had the data (**see table 1,2**) where the times are measured in years and the populations are measured in hundreds. It is common to model populations using exponentials models. 
    
    t|1|2|3|4|5
    -|--|--|--|--|--
    y|3|4|6|11|20
    
    If least squares were used to select the parameters , then we would solve 

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/leastsquare.PNG)

    t|0|2|4|6|8|10|12|14|16|18|20
    --|--|--|--|--|--|--|--|--|--|--|--
    f(t)|0|3.55|3.82|2.98|2.32|1.48|1.02|0.81|0.41|0.42|0.15

    If least squares were used to select the parameters , then we would solve 

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/new%20model.PNG)

  - **Gauss Newton**

    There are numbers of specialized methods for non linear least squares data fitting, in this module we're provide a simplest method called gauss newton. These methods compute the search direction using the the newton formula for newton methods. Suppose we have an experimental data and we want find the solution of the following nonlinear least squares problem 

    t|0|2|4|6|8|10|12|14|16|18|20
    --|--|--|--|--|--|--|--|--|--|--|--
    f(t)|0|3.55|3.82|2.98|2.32|1.48|1.02|0.81|0.41|0.42|0.15

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/new%20model.PNG)

    Suppose we want to using this module to solve nonlinear least squares problem by gauss newton method, the following ilustration is intended to solve nonlinear least squares problem with the given model above

    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::Matrix<double, 1, 11> t_i;
	    t_i << 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0;;
	    Eigen::Matrix<double, 1, 11> y_i;
	    y_i << 0.0, 3.55, 3.82, 2.98, 2.32, 1.48, 1.02, 0.81, 0.41, 0.42, 0.15;
        Eigen::Matrix<double, 1, 4> X;
	    X << -12.6218, -0.362183, 13.6189, -0.208505;
        std::cout<<"Result : "<<SK6091::functionTest::gaussNewton(X,y_i,t_i);
    }
    ```

    After compile the program : 
    > Result : -13.6218, -0.462183, 13.6189, -0.208505;

    **The Gauss Newton Method is consist of 3 parameters:**
    * *Initial Point (as a matrix type)*
    * *data population t (as matrix type)*
    * *data population f(t) as matrix type*

    The following graph give us knowledge how robust our module to estimate the nonlinear least square problem

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/gauss%20graph.jpg)

  - **Levenberg–Marquardt** 

    In this module/ library we also provide a well known algorithm for nonlinear least squares problems called Levenberg-Marquardt(LM). The LMA is used in many software applications for solving generic curve-fitting problems. However, as with many fitting algorithms, the LMA finds only a local minimum, which is not necessarily the global minimum. The LMA interpolates between the Gauss–Newton algorithm (GNA) and the method of gradient descent. The LMA is more robust than the GNA, which means that in many cases it finds a solution even if it starts very far off the final minimum. These methods compute the search direction using the damped version of newton formula for newton methods. Suppose we have an experimental data and we want find the solution of the following nonlinear least squares problem 
    
    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/fix%20damped.png)

    t|0|2|4|6|8|10|12|14|16|18|20
    --|--|--|--|--|--|--|--|--|--|--|--
    f(t)|0|3.55|3.82|2.98|2.32|1.48|1.02|0.81|0.41|0.42|0.15

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/new%20model.PNG)

    Suppose we want to using this module to solve nonlinear least squares problem by Levenberg–Marquardt method, the following ilustration is intended to solve non linear least squares problem with the given model above

    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::Matrix<double, 1, 11> t_i;
	    t_i << 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0;;
	    Eigen::Matrix<double, 1, 11> y_i;
	    y_i << 0.0, 3.55, 3.82, 2.98, 2.32, 1.48, 1.02, 0.81, 0.41, 0.42, 0.15;
        Eigen::Matrix<double, 1, 4> X;
	    X << 14.6189, -0.108505, -13.6218, -0.262183;
        std::cout << "Result : "<<SK6091::functionTest::Levenberg(X, y_i, t_i) << std::endl;

    }
    ```

    After compile the program : 
    > Result : -6.10177, -0.902389, 6.15153, -0.135121

    **The Levenberg-Marquardt Method is consist of 3 parameters:**
    * *Initial Point (as a matrix type)*
    * *data population t (as matrix type)*
    * *data population f(t) as matrix type*

    The following graph give us knowledge how robust our module to estimate the nonlinear least square probles

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/benchamark.jpg)

## 5.Calculating Derivatives
  - **Approximating the Gradient**
    
    Direct search methods are useful when the derivative of the optimization function is not available to effectively guide the search for the optimum. While direct search methods explore the parameter space in a systematic manner, they are not computationally very efficient. On the other hand, gradient methods use information from the derivatives of the optimization function to more effectively guide the search and find optimum solutions much quicker. In this module/ library we also provided a special function to aproximate a gradient of an objective function
    . The following code is intended to demonstrate how we can call the grad function to aproximate  Gradient of Rastrign function
    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::RowVector2d initalPoint;
	    initalPoint << 0.5,-0.5;
	    SK6091::functionTest::grad(initalPoint);
        return 0;
    }
    ```

    After compile the program : 
    > The Gradient of Rastrign Function :  1.000 -1.000

   
    The following benchmark results have been generated using some numerical experiments to aproximate a gradient of an objective function.

    
    N |Objective Function                 | Aproximation |Initial Point  
    --|---------------------------|--------------------------|--------------
    1|Griewank       |  0.449,   -0.215|0.5,  -0.5
    2|Rosenbrock|  149.0,   -150.0|0.5,  -0.5
    3|Ackley|  1.80967,   -1.809|0.5,  -0.5
    4|Booth_s|  -33.00,   -39.00|0.5,  -0.5
    5|Michalewicz|  -3.9e-2,   -3.8e-|0.5,  -0.5
    6|Spring_sys|  -28.34,   -180.6|0.5,  -0.5
    7|Rastrigin|  1.00,   -1.00|0.5,  -0.5

  - **Approximating the Hessian**

    The Hessian matrix or Hessian is a square matrix of second-order partial derivatives of a scalar-valued function, or scalar field. It describes the local curvature of a function of many variables. The following code is intended to demonstrate how we can call the hessian function to aproximate  hesian matrix of Rastrign function
    
    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
	    Eigen::RowVector2d initalPoint;
	    initalPoint << 0.5,-0.5;
	    SK6091::functionTest::hessian(initalPoint);
	    return 0;
    }
    ```

    After compile the program : 
    > The Hessian of Rastrign Function :  -392.783|-2.15611e+06|-2.08095e+06|-392.783

    The following code is intended to demonstrate how we can measure a computational time to aproximate hessian matrix of Rastrign function
    
    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
	    auto start = std::chrono::steady_clock::now();
	    Eigen::RowVector2d initalPoint;
	    initalPoint << 0.5,-0.5;
	    SK6091::functionTest::hessian(initalPoint);
	    auto end = std::chrono::steady_clock::now();
	    std::chrono::duration<double> sec = end - start;
	    std::cout << "Elapsed time \t:" << sec.count() << std::endl;
	    return 0;
    }
    ```

    After compile the program : 
    > Elapsed time    :0.0021832

    The following benchmark results have been generated using some numerical experiments to aproximate a hessian of an objective function.

    N |Objective Function       | f_xx | f_xy |f_yx| f_yy    
    --|-------------------------|-----------------|-----------------|---------------|--------------
    1|Griewank|0.823303|0.117375|0.117375|0.412151
    2|Rosenbrock|502|-200|-200|200
    3|Ackley|-5.63295|1.99064|1.99064|-5.63295
    4|Booth_s|10|8|8|10
    5|Michalewicz|-3.16736e-19|0|0|3.0761e-13
    6|Spring_sys|136.121|107.269|107.269|227.194
    7|Rastrigin|-392.783|-2.15611e+06|-2.08095e+06|-392.783

## 6.Heuristic-Optimization
  - **Spiral Dynamic Optimization**
    
    The spiral optimization (SPO) algorithm is a metaheuristic inspired by spiral phenomena in nature, the first SPO algorithm was proposed for two-dimensional unconstrained optimization based on two-dimensional spiral models. In this module we're also provide a brief overview how to solve constrained optimization problem. The SPO algorithm is a multipoint search algorithm that has no objective function gradient, which uses multiple spiral models that can be described as deterministic dynamical systems.The SPO algorithm for a minimization (R-2) problem under the maximum iteration (termination criterion) is as follows :  
    
    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/spiral%20(3).PNG)
    
    The following code is intended to demonstrate how we can call the Spiral Dynamic Optimization method  to minimize Griewank function through this module *(Note : namespace for this module is **SK6091** and folowed by it's class)*  
       
    ```c++
    #include "multiDimension/SpiralImp.hpp"
    int main()
    {
        std::cout<< "The Minimum Point occurs at x,y \t:\n"<<
            SK6091::Spiral::t_solve(500, 2)<<std::endl;
	    return 0;
    }
    ```

    After compile the program : 
    > The Minimum Point occurs at x,y :  1.0000 1.000

    **The Following code is Spiral Dynamic Optimization abstraction**

    ```c++
    #pragma once
    #ifndef SPIRAL_HEADER
    #define SPIRAL_HEADER
    #include <cmath>
    #include <math.h>
    #include <eigen/dense>
    #include "spiralImp.h"
    namespace SK6091{
        class Spiral {
            private:
            public:
                static Eigen::RowVector2d t_center(Eigen::Matrix<double, 50, 2>, int);
		        static Eigen::RowVector2d t_move(Eigen::RowVector2d, Eigen::Matrix<double, 2, 2>, 
                        Eigen::Matrix<double, 2, 2>, Eigen::RowVector2d);
		        static Eigen::RowVector2d t_solve(unsigned int,int);
        };
    }
    #endif
    ``` 
    **The Spiral Dynamic Optimization Method is consist of 2 parameters:**
    * *Maximum Iteration (as an unsigned int)*
    * *Dimension of Spiral Optimization Algorithm (as an unsigned int)*

    For constrained optimization problem, our module uses a penalty function to deal with some constrains. Suppose we have the nonlinear constrained optimization problem which has formulated as:

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/nonlinear%20formulation.PNG)
    
    Having transform the constrained optimization problem to unconstrained optimization problem, the general spiral optimization  algorithm for solving optimization  is as follow: 
    
    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/rn.PNG)

    The following code is intended to demonstrate how we can call the Spiral Dynamic Optimization method  for solving constrained optimization problems through this module *(Note : namespace for this module is **SK6091** and folowed by it's class)*  
       
    ```c++
    #include "multiDimension/Spiralimp.hpp"
    inline  double SK6091::Spiral::constrainedFun(Eigen::RowVector2d x, int n){
        double fx,c1;
        fx = std::pow(x[0], 2.0) + std::pow(x[1] - 1.0, 2.0);
	    c1 = x[1] - std::pow(x[0], 2.0);
	    double sumConstrained = std::pow(c1, 2.0); 
        //handling constrains, problem 5 from our paper nonlinear constrained optimization
	    fx = fx + sumConstrained;
        return fx;
    }
    int main()
    {
        std::cout<< "The Minimum Point occurs at  \t:\n"<<
            SK6091::Spiral::t_solve(500, 2)<<std::endl;
	    return 0;
    }
    ```

    After compile the program : 
    > The Minimum Point occurs at x,y :  0.001167, 0.50000

    In our paper (**Finding Solution Of Nonlinear Constrained Optimization Using Spiral Optimization Algorithm**), we try to  evaluate the effectiveness and the efficiency of the proposed method, with a given the set of nonlinear constrained optimization problems from various benchmarks problems have been tested. We summarize the result as follow (*section 4.2* **Finding Solution Of Nonlinear Constrained Optimization Using Spiral Optimization Algorithm**)

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/table%201.PNG)

## 7.Optimization-Nonlinear Equations
  - **Introduction** 

    A system of nonlinear equations is a system of two or more equations in two or more variables containing at least one equation that is not linear. In this module we're provide some robust algorithm like newton, spiral dynamic optimization, broyden method and homotopy for solving n system of nonlinear equation. The following nonlinear equation is intended to test each algorithm   

    ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/assert/fig1.PNG)

    **The Following code is objective function abstraction**
    ```c++
    #ifndef SPECIAL_FUNCTION
    #define SPECIAL_FUNCTION
    #include <cmath>
    #include <eigen/dense>
    #include <fstream>
    namespace SK6091 {
	    class functionTest
	    {
	        public:
		        inline static double F(int, Eigen::RowVector3d);
	        private:
	    };
    }
    #include "specialFunctionImp.hpp"
    #endif
    ```

    **The objective function is consist of 2 parameters:**
    * *index of element (as an integer type)*
    * *initial value (as row vector type)*

  - **Newton's Method**

    In this module we're provide a well known method called Newton method. Newton’s method is an algorithm for finding the roots of diferentiable functions, that uses iterated local linearization of a function to approximate its roots. Newton’s method also extends to systems of n differentiable functions in n variables.
    Suppose we want to using this module to solve nonlinear system of equations problem by newton method, the following ilustration is intended to solve nonlinear system of equation problem (**see figure 2.**)

    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::RowVector3d X;
	    X<< 0.1, 0.1, -0.1;
        std::cout << "Result : "<<SK6091::functionTest::nonLinearNewton(X) << std::endl;

    }
    ```

    After compile the program : 
    > Result : -6.10177, -0.902389, 6.15153

    **The Newton Method is consist of 1 parameter:**
    * *Initial Point (as a row vector type)*
      
  - **Steepest Descent**
  
    The classical steepest descent method is one of the oldest methods for the minimization of a general nonlinear function. The steepest descent method, also known as the gradient descent method, was first proposed by Cauchy in 1847 [1]. In the original paper, Cauchy proposed the use of the gradient as a way of solving a nonlinear equation. 
    Suppose we want to using this module to solve nonlinear system of equations problem by steepest descent method, the following ilustration is intended to solve nonlinear system of equation problem (**see figure 2.**)

    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::RowVector3d X;
	    X<< 0.1, 0.1, -0.1;
        std::cout << "Result : "<<SK6091::functionTest::nonLinearSteepestDescent(X) << std::endl;

    }
    ```

    After compile the program : 
    > Result : -6.10177, -0.902389, 6.15153

    **The Steepest Descent Method is consist of 1 parameter:**
    * *Initial Point (as a row vector type)*


  - **Broyden's Method**

    One of the drawbacks of using Newton’s Method to solve a system of nonlinear equations **g(x) = 0** is the computational expense that must be incurred during each iteration to evaluate the partial derivatives of g at **x(k)**, and then solve a system of linear equations involving the resulting Jacobian matrix. The algorithm does not facilitate the re-use of data from previous iterations, and in some cases evaluation of the partial derivatives can be unnecessarily costly. In this module we're also provide a method called Broyden's method that solve computational expense of using newton method. 
    Suppose we want to using this module to solve nonlinear system of equations problem by Broyden method, the following ilustration is intended to solve nonlinear system of equation problem (**see figure 2.**)

    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::RowVector3d X;
	    X<< 0.1, 0.1, -0.1;
        std::cout << "Result : "<<SK6091::functionTest::homotopyt(X) << std::endl;

    }
    ```

    After compile the program : 
    > Result : -6.10177, -0.902389, 6.15153

    **The Newton Method is consist of 1 parameter:**
    * *Initial Point (as a row vector type)*

  - **Homotopy** 

    A numerical continuation is an algorithm which takes as input a system of parametrized nonlinear equations and an initial solution and produces a set of points on the solution component. Suppose we want to using this module to solve nonlinear system of equations problem by homotopy/continuation method, the following ilustration is intended to solve nonlinear system of equation problem (**see figure 2.**)
 
    ```c++
    #include "specialFunctionImp.hpp"
    int main()
    {
        Eigen::RowVector3d X;
	    X<< 0.1, 0.1;
        std::cout << "Result : "<<SK6091::functionTest::homotopyt(X) << std::endl;

    }
    ```

    After compile the program : 
    > Result : -6.10177, -0.902389

    **The Homotopy/ continuation Method is consist of 1 parameter:**
    * *Initial Point (as a row vector type)*
