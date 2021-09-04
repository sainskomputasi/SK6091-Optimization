#pragma once
#ifndef LINES_ALGORITHM
#define LINES_ALGORITHM
#include <fstream>
namespace MA5171 {
	class Optimization {
	public:
		inline double f(double);
		inline double gradient(double);
		inline std::ofstream& write(std::ofstream&, std::string);
		inline double zoom(double, double);
		double LineSearch(double, double, double);
	private:
	};
}
#endif