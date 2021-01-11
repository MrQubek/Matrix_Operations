
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

template <>
bool MyAlgebra::CMatrix<int>::comparisionOperation(const  MyAlgebra::CMatrix<int>& other) const {
	
	if (this->rowCount != other.rowCount && this->columnCount != other.columnCount) {
		return false;
	}
	
	for (int i = 0, j = 0; i < this->rowCount; i++) {
		for (j = 0; j < this->columnCount; j++) {
			if (rowPtr[i][j] != other.rowPtr[i][j])
				return false;
		}
	}
	return true;
}

template <>
bool MyAlgebra::CMatrix<float>::comparisionOperation(const  MyAlgebra::CMatrix<float>& other) const {

	if (this->rowCount != other.rowCount && this->columnCount != other.columnCount) {
		return false;
	}

	for (int i = 0, j = 0; i < this->rowCount; i++) {
		for (j = 0; j < this->columnCount; j++) {
			if (abs(this->rowPtr[i][j] - other.rowPtr[i][j]) > ALG_PRECISION){
				return false;
			}
		}
	}
	return true;
}

template <>
bool MyAlgebra::CMatrix<double>::comparisionOperation(const  MyAlgebra::CMatrix<double>& other) const {

	if (this->rowCount != other.rowCount && this->columnCount != other.columnCount) {
		return false;
	}

	for (int i = 0, j = 0; i < this->rowCount; i++) {
		for (j = 0; j < this->columnCount; j++) {
			if (abs(rowPtr[i][j] - other.rowPtr[i][j]) > ALG_PRECISION)
				return false;
		}
	}
	return true;
}
