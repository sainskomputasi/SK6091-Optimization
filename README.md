# SK6091-Optimization
  **SK6091-Optim** is a lightweight C++ library of numerical optimization methods, it is written by Michael sipayung. The goals are to complete an independent research project and to build an Optimization library for C++ which can provide solutions for a single variable and multi-variable optimization. The users should be able to select the algorithm to be used and also access any class abstraction. Available Algorithms: golden section search, dichotomous Search, Fibonacci Search, Steepest-Descent Method, Newton Method, Gauss-Newton Method, quasi newton method, Particle Swarm Optimization, Simulated annealing, Spiral Dynamic Optimization, Differential Evolution, Evolutionary algorithms. I also provide a benchmark for each algorithm at the end of each section that consists of the optimization algorithm. The following graph is intended to be used as a sample to compile. 
  ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/my%20head.png)
  ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/ezgif.com-gif-maker%20(6).gif)
  

## 1.The optimization Problem
  - **Introduction**
  - **The Basic Optimization Problem** 
  - **How to compile** 
## 2.One Dimensional Optimization
  - **Introduction**
  - **Newton Method**
    - Newton’s method is an open method similar to Newton-Raphson because it does not require initial guesses that bracket the optimum [a,b]. In the newton method in this library, we need not a first derivative and second derivative to calculate, because we already provide an approximation to them. See The finite-difference derivative approximations for calculating derivatives in the calculating derivative section.The following code is intended to demonstrate how we call the golden section method through this module. *(Note : namespace for this module is **SK6091** and folowed by it's class)*.   
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

        **The Following code is golden section abstraction**
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

  - **Fibonacci Search**
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
	        std::cout << "maximum occurs at x \t: " << testSecant.secant(objectiveFunc,0.0,4.0,0.00001, 100);
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
		        double secant(double(double), double,double, double, int);
	        private:
	        };
        }
        #include "oneDimensional.hpp"
        #endif
        ```

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
