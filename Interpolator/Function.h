#pragma once
#include <memory>
#include <vector>
#include <tuple>
#include <list>

typedef std::vector<double> Polynomial;
typedef std::tuple<std::shared_ptr<Polynomial>, double, double> Spline;

class Function
{
	public:
		Function();
		Function(std::shared_ptr<Polynomial> p, double start = 0, double end = 0);
		Function& addSpline(std::shared_ptr<Polynomial> p, double start = 0, double end = 0);
		double operator()(double x) const;
		~Function();

	private:
		std::list<Spline> splines;
};

