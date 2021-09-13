# SK6091-Optimization
  SK6091-Optim adalah sebuah module yang berguna untuk menyelesaikan permasalahan-permasalahan optimasi. Algoritma yang dikembangkan saat ini masih memiliki keterbatasan jumlah algoritma, module ini dikembangkan dengan bahasa c++ yang hanya berisikan header/ abstract dan implementasi dari algoritma yang dapat dijalankan pada compiler yang mendukung c++ versi c-11 . Module ini juga berguna untuk melengkapi riset perkuliahan yang diberikan secara mandiri. Berikut merupakan contoh-contoh fungsi yang dapat digunakan dalam menyelesaikan permasalahan optimasi Numerik. 
  ![This is an image](https://upload.wikimedia.org/wikipedia/commons/7/72/Max_paraboloid.svg)
  
  _Graph dari fungsi  f(x, y) = −(x² + y²) + 4 _

## 1.The optimization Problem
  - Introduction
  - The Basic Optimization Problem 
  - How to compile 
## 2.One Dimensional Optimization
  - Introduction
  - Dichotomous Search
  - Fibonacci Search
  - Golden-Section Search
    - The golden-section search is a simple, general-purpose, single-variable search technique. It is similar in spirit to the bisection approach for locating root. The following code is intended to demonstrate how we call the golden section method through this module. (Note : namespace for this module is SK6091 and folowed by it's class) 
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
        The Following code is golden section abstraction
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
    
  - Quadratic Interpolation Method
  - Cubic Interpolation 
## 3.Multi-Dimensional Gradient Method
  - Introduction
  - Steepest-Descent Method
  - Newton Method
  - Gauss Newton Method
## 4.Quasi-Newton Method 
  - Introduction 
  - Generating Matrix S
  - Rank One Method 
  - Hoshino Method
  - Davidon-Fletcher-Powell Method
## 5.Calculating Derivatives
  - Finite-Difference Derivative Approximations
  - Approximating the Gradient
  - Approximating a Sparse Jacobian
  - Approximating the Hessian
  - Approximating a Sparse Hessian
## 6.Heuristic-Optimization
  - Particle Swarm Optimization
  - Simulated annealing
  - Spiral Dynamic Optimization
  - Differential evolution
  - Evolutionary algorithms
## 7.Nonlinear Equations
  - Newton's Method
  - Inexact Newton's Method
  - Broyden's Method
  - Tensor Method 
