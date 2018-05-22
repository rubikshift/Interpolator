#include <iostream>
#include <fstream>
#include <locale>
#include "Interpolator.h"

struct comma_separator : std::numpunct<char>
{
	virtual char do_decimal_point() const override
	{
		return ',';
	}
};

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "Za mala liczba parametrow :(" << std::endl;
		return -1;
	}

	Interpolator interpolator;
	interpolator.loadData(argv[1]).choseInterpolationParts(5);
	
	{
		auto lagrangeFunction = interpolator.lagrange();
		auto dataSet = interpolator.getDataSet();
		std::ofstream output("langrange.csv", std::ios::out);
		output.imbue(std::locale(output.getloc(), new comma_separator));
		output << "x; y" << std::endl;
		for(auto const& p : dataSet)
			output << p.first << "; " << lagrangeFunction(p.first) << std::endl;
		output.close();
	}

	{
		auto splineFunction = interpolator.spline();
		auto dataSet = interpolator.getDataSet();
		std::ofstream output("spline.csv", std::ios::out);
		output.imbue(std::locale(output.getloc(), new comma_separator));
		output << "x; y" << std::endl;
		for (auto const& p : dataSet)
			output << p.first << "; " << splineFunction(p.first) << std::endl;
		output.close();
	}
	return 0;
}