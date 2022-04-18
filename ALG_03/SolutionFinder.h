#pragma once
#include <vector>
class SolutionFinder
{
private:


	// ��������� ������������ ��������.
	double findScalar(const std::vector<double>& x1, const std::vector<double>& x2);

	// ����� ������ �������.
	void findX(const std::vector<std::vector <double>>& mAB, std::vector<double>& x);

	// ���������� ������� � ������������ ����.
	void triangularMatrix(std::vector<std::vector <double>>& mAB);


	// ���������� �������� �������.
	void inverseMatrix(const std::vector<std::vector <double>>& mA, 
		std::vector<std::vector <double>>& inMatrix);

	// A = A + E*value
	void findHelpMatrix(const std::vector<std::vector <double>>& mA,
		std::vector<std::vector <double>>& mB, double value);

public:

	double findMaxAbsEigenvalue(const std::vector<std::vector <double>>& mA,
		const double& e, std::vector<double> x, double givenEigen);

	
	double findMaxEigenvalue(const std::vector<std::vector <double>>& mA, double e,
		const double& maxAbsEigenvalue, std::vector<double> x, double givenEigen);


	double findMinEigenvalue(const std::vector<std::vector <double>>& mA, double e,
		double maxAbsEigenvalue, std::vector<double> x, double givenEigen);
		

	double findClosestEigenvalue(std::vector<std::vector <double>>& mA, double e,
		std::vector<double> x, double givenEigenvalue);







	



};

