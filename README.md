# SK6091-Optimization
  **SK6091-Optim** is a lightweight C++ library of numerical optimization methods, it is written by Michael sipayung. The goals are to complete an independent research project and to build an Optimization library for C++ which can provide solutions for a single variable and multi-variable optimization. The users should be able to select the algorithm to be used and also access any class abstraction. Available Algorithms: golden section search, dichotomous Search, Fibonacci Search, Steepest-Descent Method, Newton Method, Gauss-Newton Method, quasi newton method, Particle Swarm Optimization, Simulated annealing, Spiral Dynamic Optimization, Differential Evolution, Evolutionary algorithms. I also provide a benchmark for each algorithm at the end of each section that consists of the optimization algorithm. The following graph is intended to be used as a sample to compile. 
  ![This is an image](https://github.com/sainskomputasi/SK6091-Optimization/blob/master/NumericalOptimization/my%20head.png)
  

## 1.The optimization Problem
  - **Introduction**
  - **The Basic Optimization Problem** 
  - **How to compile** 
## 2.One Dimensional Optimization
  - **Introduction**
  - **Dichotomous Search**
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
	        std::cout<<"Minimum Value\t: "<< testGolden.goldenSec(objectiveFunc, 0.0, 4.0, 100, 0.00001);
        }
        ```
        After compile the program : 
        > Minimum Value : 1.42755

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
  - Quadratic Interpolation Method
  - Cubic Interpolation 
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
