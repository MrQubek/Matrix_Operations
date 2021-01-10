
#include "CMatrix.h"

template <>
void MyAlgebra::CMatrix<int>::populateMatrixWithRandomNumbers() {
	static int seed = clock() * time(NULL);
	std::default_random_engine rndNrGenerator;
	rndNrGenerator.seed(seed++);
	std::uniform_int_distribution<int> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);

	for (int i = 0, j; i < rowCount; i++) {
		for (j = 0; j < columnCount; j++) {
			rowPtr[i][j] = distribution(rndNrGenerator);
		}
	}
}

template <>
void MyAlgebra::CMatrix<float>::populateMatrixWithRandomNumbers() {
	static int seed = clock() * time(NULL);
	std::default_random_engine rndNrGenerator;
	rndNrGenerator.seed(seed++);
	std::uniform_real_distribution<float> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);

	for (int i = 0, j; i < rowCount; i++) {
		for (j = 0; j < columnCount; j++) {
			rowPtr[i][j] = distribution(rndNrGenerator);
		}
	}
}

template <>
void MyAlgebra::CMatrix<double>::populateMatrixWithRandomNumbers() {
	static int seed = clock() * time(NULL);
	std::default_random_engine rndNrGenerator;
	rndNrGenerator.seed(seed++);
	std::uniform_real_distribution<double> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);

	for (int i = 0, j = 0; i < rowCount; i++) {
		for (j = 0; j < columnCount; j++) {
			rowPtr[i][j] = distribution(rndNrGenerator);
		}
	}
}