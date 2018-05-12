#pragma once
#include <memory>
#include <vector>
#include <tuple>
#include "Function.h"

typedef std::pair<double, double> Point;
typedef std::vector<Point> DataVector;
typedef std::vector<double> Polynomial;

class Interpolator
{
	public:
		Interpolator();
		Interpolator& loadData(std::string dataSetFileName);
		Interpolator& choseInterpolationPoints(std::size_t points);
		Interpolator& randomlyChooseInterpolationPoints(std::size_t points);
		Function lagrange();
		Function spline();
		const DataVector& getDataSet() const;
		~Interpolator();

	private:
		DataVector dataSet;
		DataVector interpolationSet;
		std::size_t points;

		static std::shared_ptr<Polynomial> multiplyPolynomials(std::shared_ptr<Polynomial> fi1, std::shared_ptr<Polynomial> fi2);
		static void addPolynomials(std::shared_ptr<Polynomial> dst, const std::shared_ptr<Polynomial>& src);

};

