#include <iostream>
#include <Eigen/Dense>
#include <dlib/matrix.h>
#include <dlib/optimization.h>

using namespace Eigen;

struct calc
{
	double operator() (const dlib::matrix<double, 0, 1> x) const
	{
		return x(0)*x(0)*x(0);
	}
};

struct derivate
{
	dlib::matrix<double, 0, 1> operator() (const dlib::matrix<double, 0, 1> x) const
	{
		dlib::matrix<double, 0, 1> de;
		de.set_size(x.size());
		for (int i = 0; i < x.size(); i++)
			de(i) = 3*x(i)*x(i);
		return de;
	}
};



int main()
{
	calc _calc;
	derivate _derivate;
	double E;
	dlib::matrix<double,0,1> x;
	x.set_size(1);
	x(0) = 1;

	E = dlib::find_min(dlib::bfgs_search_strategy(),
		dlib::objective_delta_stop_strategy(1e-6f, 1000),
		_calc,
		_derivate,
		x,
		-10.0);

	std::cout << E<<" "<<x;

}
