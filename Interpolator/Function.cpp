#include <cmath>
#include "Function.h"

Function::Function(std::shared_ptr<Polynomial> p, double start, double end)
{
	splines.emplace_back(p, start, end);
}

Function& Function::addSpline(std::shared_ptr<Polynomial> p, double start, double end)
{
	splines.emplace_back(p, start, end);
	return *this;
}

double Function::operator()(double x) const
{
	for (auto const& f : splines)
	{
		auto polynomial = std::get<0>(f);
		auto start = std::get<1>(f);
		auto end = std::get<2>(f);

		if (start != end && x >= start && x <= end)
			continue;

		double y = 0;
		for (Polynomial::size_type i = 0; i < polynomial->size(); i++)
			y += polynomial->at(i) * pow(x, polynomial->size() - i - 1);
		return y;
	}
	return NAN;
}

Function::~Function()
{
}
