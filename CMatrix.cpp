
#include "CMatrix.h"

template <>
void MyAlgebra::CMatrix<int>::populateMatrixWithRandomNumbers() {
	for (int i = 0, j = 0; i < rowCount; i++) {
		for (j = 0; j < rowCount; j++) {
			rowPtr[i][j] = rand() % (RANDOM_MAX_VALUE - RANDOM_MIN_VALUE) + RANDOM_MIN_VALUE;
		}
	}
}

template <>
void MyAlgebra::CMatrix<float>::populateMatrixWithRandomNumbers() {

	std::default_random_engine rndNrGenerator;
	std::uniform_real_distribution<float> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);

	for (int i = 0, j = 0; i < rowCount; i++) {
		for (j = 0; j < rowCount; j++) {
			rowPtr[i][j] = distribution(rndNrGenerator);
		}
	}
}

template <>
void MyAlgebra::CMatrix<double>::populateMatrixWithRandomNumbers() {

	std::default_random_engine rndNrGenerator;
	std::uniform_real_distribution<double> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);

	for (int i = 0, j = 0; i < rowCount; i++) {
		for (j = 0; j < rowCount; j++) {
			rowPtr[i][j] = distribution(rndNrGenerator);
		}
	}
}