#include <fstream>
#include <locale>
#include <algorithm>
#include <cmath>
#include <random>
#include "Interpolator.h"

#include "matrix.h"
#include "solvingmethods.h"

Interpolator::Interpolator()
{
}

Interpolator& Interpolator::loadData(std::string dataSetFileName)
{
	std::ifstream input(dataSetFileName, std::ios::in);

	double x, y;

	while (input >> x >> y)
		dataSet.emplace_back(x, y);

	input.close();
	return *this;
}

Interpolator& Interpolator::choseInterpolationParts(std::size_t parts)
{
	if (parts == 0)
		throw std::invalid_argument("Ilosc przedzialow nie moze byc mniejsza rowna 0");
	else if (parts + 1 > dataSet.size())
		throw std::invalid_argument("Zbyt duza liczba przedzialow");

	interpolationSet.clear();

	points = parts + 1;
	std::size_t dx = dataSet.size() / points;
	if (dataSet.size() != dx * points)
		dx++;

	for (DataVector::size_type i = 0; i < dataSet.size(); i += dx)
		interpolationSet.emplace_back(dataSet[i]);
	if (dataSet.size() != dx * points)
	{
		interpolationSet.pop_back();
		interpolationSet.emplace_back(dataSet.back());
	}
	return *this;
}

Interpolator& Interpolator::randomlyChoseInterpolationParts(std::size_t parts)
{
	if (parts == 0)
		throw std::invalid_argument("Ilosc przedzialow nie moze byc mniejsza rowna 0");
	else if (parts + 1 > dataSet.size())
		throw std::invalid_argument("Zbyt duza liczba przedzialow");

	interpolationSet.clear();

	this->points = parts + 1;

	auto xs = new int[points];

	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, dataSet.size() - 1);
	
	int random = 0;
	bool ok = false;
	for (std::size_t i = 0; i < points; i++)
	{
		ok = false;
		while (!ok & i > 0)
		{
			random = distribution(generator);
			for (std::size_t j = 0; j < i; j++)
			{
				ok = false;
				if (xs[j] == random)
					break;
				ok = true;
			}
		}
		xs[i] = random;
	}

	std::sort(xs, xs + points);

	for (std::size_t i = 0; i < points; i++)
		interpolationSet.emplace_back(dataSet[xs[i]]);

	delete xs;
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
	const std::size_t equations = 4 * (points - 1);
	auto matrix = std::make_shared<Matrix>(equations, equations, 0);
	auto left = std::make_shared<Matrix>(equations, 1); // [d0 c0 b0 a0 d1 c1 b1 a1 .. di ci bi ai]^T

	//hi = x(i+1) - xi
	auto h = std::make_unique<Matrix>(points - 1, 1);
	for (std::size_t i = 1; i < points; i++)
		h->set(i - 1, interpolationSet[i].first - interpolationSet[i - 1].first);

	//c0 = 0
	matrix->set(1, 1, 1);
	left->set(1, 0);

	//ai = f(xi)
	std::size_t helpIter = 0;
	for (std::size_t i = 3; i < equations; i += 4)
	{
		matrix->set(i, i, 1);
		left->set(i, interpolationSet[helpIter].second);
		helpIter++;
	}

	//di * 6hi + 2ci - 2c(i+1) = 0
	helpIter = 0;
	for (std::size_t i = 0; i < equations; i += 4)
	{
		matrix->set(i, i, 6 * h->at(helpIter));	//di * 6hi 
		matrix->set(i, i + 1, 2);				//2ci
		if (i + 5 < equations)
			matrix->set(i, i + 5, -2);			//-2c(i+1)
		left->set(i, 0);
		helpIter++;
	}

	//di * hi^3 + ci * hi^2 + bi * hi + ai = f(x(i+1))
	helpIter = 0;
	for (std::size_t i = 2; i < equations; i += 4)
	{
		matrix->set(i, 0 + 4 * helpIter, pow(h->at(helpIter), 3));	//di * hi^3
		matrix->set(i, 1 + 4 * helpIter, pow(h->at(helpIter), 2));	//ci * hi^2
		matrix->set(i, 2 + 4 * helpIter, h->at(helpIter));			//bi * hi
		matrix->set(i, 3 + 4 * helpIter, 1);						//ai
		left->set(i, interpolationSet[helpIter + 1].second);
		if (i == 2)
			i--; //only for S0(x) this equation is 3rd, for Si(x) it is 2nd
		helpIter++;
	}

	//di * 3hi^2 + ci * 2hi + bi - b(i+1) = 0; i > 0
	helpIter = 0;
	for (std::size_t i = 6; i < equations; i += 4)
	{
		matrix->set(i, i - 6, 3*pow(h->at(helpIter), 2));	//di * 3hi^2 
		matrix->set(i, i - 6 + 1, 2*h->at(helpIter));		//ci * 2hi
		matrix->set(i, i - 6 + 2, 1);						//bi
		matrix->set(i, i - 6 + 6, -1);						//- b(i+1)
		left->set(i, 0);
	}

	auto x = std::make_shared<Matrix>(matrix->rows, left->cols, 1);
	LUdecomposition(matrix, left, x);

	Function f;
	auto fi = std::make_shared<Polynomial>();
	helpIter = 0;
	for (std::size_t i = 0; i < equations; i++)
	{
		fi->emplace_back(x->at(i));

		if (i % 4 == 3)
		{
			f.addSpline(move(fi), interpolationSet[helpIter].first, interpolationSet[helpIter + 1].first);
			fi = std::make_shared<Polynomial>();
			helpIter++;
		}
	}
	return f;

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
