#include <fstream>
#include <locale>
#include <algorithm>
#include "Interpolator.h"

Interpolator::Interpolator()
{
}

Interpolator& Interpolator::loadData(std::string dataSetFileName)
{
	std::ifstream input(dataSetFileName, std::ios::in);

	double x, y;

	while (!input.eof())
	{
		input >> x;
		input >> y;
		dataSet.emplace_back(x, y);
	}

	input.close();
	return *this;
}

Interpolator& Interpolator::choseInterpolationPoints(std::size_t points)
{
	if (points < 2)
		throw std::invalid_argument("Ilosc punktow musi byc wieksza od 1");

	this->points = points;
	interpolationSet.reserve(points);
	std::size_t dx = dataSet.size() / points;

	for (DataVector::size_type i = 0; i < dataSet.size(); i += dx)
		interpolationSet.emplace_back(dataSet[i]);

	return *this;
}

Interpolator& Interpolator::randomlyChooseInterpolationPoints(std::size_t points)
{
	if (points < 2)
		throw std::invalid_argument("Ilosc punktow musi byc wieksza od 1");

	auto xs = std::make_unique<int*>(new int[points]);

	this->points = points;
	return *this;
}

Function Interpolator::lagrange()
{
	auto fi = std::make_shared<Polynomial>();
	for (std::size_t i = 0; i < points; i++)
		fi->emplace_back(0);

	double d;
	for (std::size_t i = 0; i < points; i++)
	{
		d = 1;
		auto fi1 = std::make_shared<Polynomial>();
		fi1->emplace_back(1);

		for (std::size_t j = 0; j < points; j++)
		{
			if (i == j)
				continue;
			auto fi2 = std::make_shared<Polynomial>();
			fi2->emplace_back(1);
			fi2->emplace_back(-1 * interpolationSet[j].first);

			fi1 = multiplyPolynomials(move(fi1), move(fi2));
			d *= interpolationSet[i].first - interpolationSet[j].first;
		}

		for (auto& a : *fi1)
		{
			a /= d;
			a *= interpolationSet[i].second;
		}

		addPolynomials(fi, fi1);		
	}

	return Function(move(fi));
}

Function Interpolator::spline()
{
	auto fi = std::make_shared<Polynomial>(points);
	for (std::size_t i = 0; i < points; i++)
		fi->emplace_back(0);
	return Function(move(fi));
}

const DataVector& Interpolator::getDataSet() const
{
	return dataSet;
}

Interpolator::~Interpolator()
{
}

std::shared_ptr<Polynomial> Interpolator::multiplyPolynomials(std::shared_ptr<Polynomial> fi1, std::shared_ptr<Polynomial> fi2)
{
	auto fi = std::make_shared<Polynomial>();
	for (std::size_t i = 0; i < fi1->size() + fi2->size() - 1; i++)
		fi->emplace_back(0);

	for (Polynomial::size_type i = 0; i < fi1->size(); i++)
		for (Polynomial::size_type j = 0; j < fi2->size(); j++)
		{
			fi->at(i + j) += fi1->at(i) * fi2->at(j);
		}

	return fi;
}

void Interpolator::addPolynomials(std::shared_ptr<Polynomial> dst, const std::shared_ptr<Polynomial>& src)
{
	if (dst->size() < src->size())
		throw std::invalid_argument("Wielomian src nie moze byc wiekszego stopinia niz dst");

	Polynomial::size_type i = dst->size() - 1;
	Polynomial::size_type j = src->size() - 1;
	while (j != (std::vector<int>::size_type) - 1 && i != (std::vector<int>::size_type) - 1)
	{
		dst->at(i) += src->at(j);
		i--;
		j--;
	}
}
