#include "Matrix.h"

Matrix::Matrix(std::vector<std::vector<double>> t_matrix, double t_e, double t_givenEigenValue,
	std::vector<double> t_x) : matrix(t_matrix), e(t_e), givenEigenValue(t_givenEigenValue),
	x(t_x) {}

const std::vector<std::vector<double>>& Matrix::getMatix() const
{
	return matrix;
}

double Matrix::getE() const { return e;}

double Matrix::getGivenEigenValue() const {	return givenEigenValue;}

const std::vector<double>& Matrix::getX() const { return x; }
