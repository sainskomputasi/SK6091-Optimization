# CMakeList.txt : CMake project for NumericalOptimization, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
 add_executable (NumericalOptimization "NumericalOptimization.cpp" "NumericalOptimization.h" "LIneSearchAlgorithm.h" "oneDimen/oneDimensional_Abstract.hpp" "specialFunction.hpp" "specialFunction.hpp" "multiDimension/multiDimen.hpp")
 find_package(Eigen3 CONFIG REQUIRED)
 target_link_libraries(NumericalOptimization PRIVATE Eigen3::Eigen)
# TODO: Add tests and install targets if needed.
