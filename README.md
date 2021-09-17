# SK6091-Optimization
  **SK6091-Optim** is a lightweight C++ library of numerical optimization methods, it is written by Michael sipayung. The goals are to complete an independent research project and to build an Optimization library for C++ which can provide solutions for a single variable and multi-variable optimization. The users should be able to select the algorithm to be used and also access any class abstraction. Available Algorithms: golden section search, dichotomous Search, Fibonacci Search, Steepest-Descent Method, Newton Method, Gauss-Newton Method, quasi newton method, Particle Swarm Optimization, Simulated annealing, Spiral Dynamic Optimization, Differential Evolution, Evolutionary algorithms. I also provide a benchmark for each algorithm at the end of each section that consists of the optimization algorithm. The following graph is intended to be used as a sample to compile. 
  ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/my%20head.png)
  ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/ezgif.com-gif-maker%20(6).gif)
  

## 1.The optimization Problem
  - **Introduction**
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
      inline double twoD_Rastrign(const T& x) {
	    return (20+ (std::pow(x,2.0)-10*std::cos(2*M_PI*x)+
                std::pow(y,2.0)-10*std::cos(2*M_PI*y)));
      }
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
  - **Steepest-Descent Method**
  - **Newton Method**
  - **Gauss Newton Method**
## 4.Quasi-Newton Method 
  - **Introduction** 
  - **Generating Matrix S**
  - **Rank One Method** 
  - **Hoshino Method**
  - **Davidon-Fletcher-Powell Method**
## 5.Calculating Derivatives
  - **Finite-Difference Derivative Approximations**
  - **Approximating the Gradient**
  - **Approximating a Sparse Jacobian**
  - **Approximating the Hessian**
  - **Approximating a Sparse Hessian**
## 6.Heuristic-Optimization
  - **Particle Swarm Optimization**
  - **Simulated annealing**
  - **Spiral Dynamic Optimization**
  - **Differential evolution**
  - **Evolutionary algorithms**
## 7.Nonlinear Equations
  - **Newton's Method**
  - **Inexact Newton's Method**
  - **Broyden's Method**
  - **Tensor Method** 
