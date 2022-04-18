#pragma once
#include <vector>
class Matrix
{
private:
	const std::vector<std::vector<double>> matrix;
	double e;
	double givenEigenValue;
	std::vector<double> x;

public:
	Matrix(std::vector<std::vector<double>> t_matrix, double t_e, double t_givenEigenValue,
		std::vector<double> t_x);

	const std::vector<std::vector<double>>& getMatix() const;

	double getE() const;

	double getGivenEigenValue() const;

	const std::vector<double>& getX() const;

};

