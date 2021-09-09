#pragma once
#ifndef LINES_ALGORITHM
#define LINES_ALGORITHM
#include <fstream>
namespace MA5171 {
	class Optimization {
	public:
		inline double f(double,double);
		inline double fak_x1(double,double);
		inline double fak_x2(double, double);
		inline std::ofstream& write(std::ofstream&, std::string);
		inline double zoom(double a_lo, double a_hi, double phi_lo, double phi_hi, double derphi_lo, double phi(double),  
			double derphi(double), double phi0, double derphi0, double c1, double c2);
		double LineSearch(double, double, double);
		double cubicInter(double , double , double , double , double , double , double );
		double quadInter(double , double , double , double , double );


	private:
	};
}
#endif