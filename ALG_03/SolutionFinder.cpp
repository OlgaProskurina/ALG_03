#include "SolutionFinder.h"




// Скалярное произведение векторов.
double SolutionFinder::findScalar(const std::vector<double>& x1, const std::vector<double>& x2) {
	double sum = 0;
	for (int i = 0; i < x1.size(); i++)
	{
		sum += x1[i] * x2[i];
	}
	return sum;
}



// Поиск максимального собственного значения по модулю.
double SolutionFinder::findMaxAbsEigenvalue(const std::vector<std::vector <double>>& mA,
	const double& e, std::vector<double> x, double givenEigen) {
	int matrixSize = mA.size();
	double norm, eigenvalue;
	int INum = 0;
	do
	{
		eigenvalue = givenEigen;
		std::vector<double> nextX(matrixSize);
		for (int i = 0; i < matrixSize; i++)
		{
			for (int j = 0; j < matrixSize; j++)
			{
				nextX[i] += mA[i][j] * x[j];
			}
		}
		givenEigen = findScalar(x, nextX);
		norm = sqrt(findScalar(nextX, nextX));
		for (int i = 0; i < matrixSize; i++)
		{
			nextX[i] = nextX[i] / norm;
		}
		givenEigen /= findScalar(x, x);
		x = nextX;
	} while (abs(eigenvalue - givenEigen) > e);
	return eigenvalue;
}


// Поиск корней системы.
void SolutionFinder::findX(const std::vector<std::vector <double>>& mAB, std::vector<double>& x) {
	int matrixSize = mAB.size();
	x[matrixSize - 1] = mAB[matrixSize - 1][matrixSize];
	for (int i = matrixSize - 2; i >= 0; i--) {
		x[i] = mAB[i][matrixSize];
		for (int j = i + 1; j < matrixSize; j++) {

			x[i] -= mAB[i][j] * x[j];
		}
	}
}


// Приведение матрицы к треугольному виду.
void SolutionFinder::triangularMatrix(std::vector<std::vector <double>>& mAB) {
	int istr, iclmn;
	double mem;
	int i, j;
	int matrixSize = mAB.size();
	for (i = 0; i < matrixSize; i++) {
		mem = abs(mAB[i][i]);
		istr = i;
		for (j = i + 1; j < matrixSize; j++)
		{
			if (abs(mAB[j][i]) > mem) {
				istr = j;
				mem = abs(mAB[j][i]);

			}
		}
		if (istr != i) {
			for (iclmn = i; iclmn < matrixSize+1; iclmn++) {
				std::swap(mAB[i][iclmn], mAB[istr][iclmn]);
			}
		}

		for (iclmn = i + 1; iclmn < matrixSize+1; iclmn++) {
			mAB[i][iclmn] = mAB[i][iclmn] / mAB[i][i];
		}
		mAB[i][i] = 1;
		for (istr = i + 1; istr < matrixSize; istr++)
		{
			if (mAB[istr][i] != 0) {
				for (iclmn = i + 1; iclmn < matrixSize+1; iclmn++) {
					mAB[istr][iclmn] = mAB[istr][iclmn] - (mAB[istr][i] * mAB[i][iclmn]);
				}
			}
			mAB[istr][i] = 0;
		}

	}
}


// Нахождение обратной матрицы.
void SolutionFinder::inverseMatrix(const std::vector<std::vector <double>>& mA, std::vector<std::vector <double>>& inMatrix) {
	int matrixSize = mA.size();
	std::vector<std::vector <double>> idmatrix(matrixSize, std::vector<double>(matrixSize));
	for (int i = 0; i < matrixSize; i++) {
		idmatrix[i][i] = 1;
	}
	std::vector<std::vector <double>> hmatrix(matrixSize, std::vector<double>(matrixSize+1));
	std::vector<double> x(matrixSize);
	for (int m = 0; m < matrixSize; m++) {
		for (int i = 0; i < matrixSize; i++) {
			for (int j = 0; j < matrixSize; j++) {
				hmatrix[i][j] = mA[i][j];
			}
			hmatrix[i][matrixSize] = idmatrix[i][m];
		}
		triangularMatrix(hmatrix);
		findX(hmatrix, x);
		for (int i = 0; i < matrixSize; i++) {
			inMatrix[i][m] = x[i];
		}
	}

}


// A = A + E*value
void SolutionFinder::findHelpMatrix(const std::vector<std::vector <double>>& mA,
	std::vector<std::vector <double>>& mB, double value) {
	mB = mA;
	for (int i = 0; i < mA.size(); i++) {
		mB[i][i] += value;
	}
}

// Максимальное собственное значение матрицы.
double SolutionFinder::findMaxEigenvalue(const std::vector<std::vector <double>>& mA, double e,
	const double& maxAbsEigenvalue, std::vector<double> x, double givenEigen) {
	int matrixSize = mA.size();
	std::vector<std::vector <double>> mB(matrixSize, std::vector<double>(matrixSize));
	double maxEigenvalue;
	findHelpMatrix(mA, mB, maxAbsEigenvalue);
	maxEigenvalue = findMaxAbsEigenvalue(mB, e, x, givenEigen);
	return maxEigenvalue - maxAbsEigenvalue;
}


//Минимальное собственное значение матрицы.
double SolutionFinder::findMinEigenvalue(const std::vector<std::vector <double>>& mA, double e,
	double maxAbsEigenvalue, std::vector<double> x, double givenEigen) {
	int matrixSize = mA.size();
	std::vector<std::vector <double>> mB(matrixSize, std::vector<double>(matrixSize));
	double minEigenvalue;
	maxAbsEigenvalue = -maxAbsEigenvalue;
	findHelpMatrix(mA, mB, maxAbsEigenvalue);
	minEigenvalue = findMaxAbsEigenvalue(mB, e, x, givenEigen);
	return minEigenvalue - maxAbsEigenvalue;
}



double SolutionFinder::findClosestEigenvalue(std::vector<std::vector <double>>& mA, double e,
	std::vector<double> x, double givenEigenvalue)
{
	int matrixSize = mA.size();
	std::vector<std::vector <double>> mB = mA;
	findHelpMatrix(mA, mB, -givenEigenvalue);
	std::vector<std::vector <double>> imB(matrixSize, std::vector<double>(matrixSize));
	inverseMatrix(mB, imB);
	double closestEigenvalue = findMaxAbsEigenvalue(imB, e, x, givenEigenvalue);
	return givenEigenvalue + 1.0f / closestEigenvalue;
}


