
#include "CMatrix.h"

template <>
void MyAlgebra::CMatrix<int>::populateMatrixWithRandomNumbers() {
	static int seed = clock() * time(NULL);
	std::default_random_engine rndNrGenerator;
	rndNrGenerator.seed(seed++);

#if LIMIT_MAX_RANDOM_VALUE !=0
	std::uniform_int_distribution<int> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);
#else
	std::uniform_int_distribution<int> distribution(0);
#endif

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

#if LIMIT_MAX_RANDOM_VALUE !=0
	std::uniform_real_distribution<float> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);
#else
	std::uniform_real_distribution<float> distribution(0);
#endif

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

#if LIMIT_MAX_RANDOM_VALUE !=0
	std::uniform_real_distribution<double> distribution(RANDOM_MIN_VALUE, RANDOM_MAX_VALUE);
#else
	std::uniform_real_distribution<double> distribution(0);
#endif

	for (int i = 0, j = 0; i < rowCount; i++) {
		for (j = 0; j < columnCount; j++) {
			rowPtr[i][j] = distribution(rndNrGenerator);
		}
	}
}

template <>
void MyAlgebra::CMatrix<int>::threadMultiplyOperation(
	CMatrix<int>& retMatrix, const CMatrix<int>& other, const int otherColumnIndexStart, const int otherColumnIndexEnd) const {

	if (otherColumnIndexStart > otherColumnIndexEnd) {
		throw(TwoWrongArguments("OP_MULTIPLY_MAT_THREAD", otherColumnIndexStart, otherColumnIndexEnd));
	}

	int i, j, k;

	const int commonDimLength = this->getColumnCount();

	int* otherColumnPtr = new int[commonDimLength];

	__m256i accVector;
	int accArray[SIMD_INT_LENGTH];

	for (int otherColumnIndex = otherColumnIndexStart; otherColumnIndex <= otherColumnIndexEnd; otherColumnIndex++) {

		for (i = 0; i < commonDimLength; i++) {
			otherColumnPtr[i] = other.rowPtr[i][otherColumnIndex];
		}

		//for evry row in matrix
		for (i = 0; i < this->rowCount; i++) {
			accVector = _mm256_setzero_si256();
			//multiply and add elems 
			for (j = 0; j < commonDimLength - commonDimLength % SIMD_INT_LENGTH; j += SIMD_FLOAT_LENGTH) {
				accVector = _mm256_add_epi32(_mm256_mullo_epi32(
					_mm256_loadu_epi32(&this->rowPtr[i][j]), _mm256_loadu_epi32(&otherColumnPtr[j])),
					accVector);
			}

			_mm256_storeu_epi32(accArray, accVector);
			for (k = 1; k < SIMD_INT_LENGTH; k++) {
				accArray[0] += accArray[k];
			}

			for (j = commonDimLength - commonDimLength % SIMD_INT_LENGTH; j < commonDimLength; j++) {
				accArray[0] += this->rowPtr[i][j] * otherColumnPtr[j];
			}
			retMatrix.rowPtr[i][otherColumnIndex] = accArray[0];
		}
	}
	delete[] otherColumnPtr;
}

template <>
void MyAlgebra::CMatrix<float>::threadMultiplyOperation(
	CMatrix<float>& retMatrix, const CMatrix<float>& other, const int otherColumnIndexStart, const int otherColumnIndexEnd) const {
	
	if (otherColumnIndexStart > otherColumnIndexEnd) {
		throw(TwoWrongArguments("OP_MULTIPLY_MAT_THREAD", otherColumnIndexStart, otherColumnIndexEnd));
	}

	int i, j, k;

	const int commonDimLength = this->getColumnCount();

	float* otherColumnPtr = new float[commonDimLength];

	__m256 accVector;
	float accArray[SIMD_FLOAT_LENGTH];

	for (int otherColumnIndex = otherColumnIndexStart; otherColumnIndex <= otherColumnIndexEnd; otherColumnIndex++) {

		for (i = 0; i < commonDimLength; i++) {
			otherColumnPtr[i] = other.rowPtr[i][otherColumnIndex];
		}

		//for evry row in matrix
		for (i = 0; i < this->rowCount; i++) {
			accVector = _mm256_setzero_ps();
			//multiply and add elems 
			for (j = 0; j < commonDimLength - commonDimLength % SIMD_FLOAT_LENGTH; j += SIMD_FLOAT_LENGTH) {
				accVector = _mm256_fmadd_ps(_mm256_loadu_ps(&this->rowPtr[i][j]), _mm256_loadu_ps(&otherColumnPtr[j]), accVector);
			}

			_mm256_storeu_ps(accArray, accVector);
			for (k = 1; k < SIMD_FLOAT_LENGTH; k++) {
				accArray[0] += accArray[k];
			}

			for (j = commonDimLength - commonDimLength % SIMD_FLOAT_LENGTH; j < commonDimLength; j++) {
				accArray[0] += this->rowPtr[i][j] * otherColumnPtr[j];
			}
			retMatrix.rowPtr[i][otherColumnIndex] = accArray[0];
		}
	}
	delete[] otherColumnPtr;
}

template <>
void MyAlgebra::CMatrix<double>::threadMultiplyOperation(
	CMatrix<double>& retMatrix, const CMatrix<double>& other, const int otherColumnIndexStart, const int otherColumnIndexEnd) const {

	if (otherColumnIndexStart > otherColumnIndexEnd) {
		throw(TwoWrongArguments("OP_MULTIPLY_MAT_THREAD", otherColumnIndexStart, otherColumnIndexEnd));
	}

	int i, j, k;
	const int commonDimLength = this->getColumnCount();

	double* otherColumnPtr = new double[commonDimLength];

	__m256d accVector;
	double accArray[SIMD_DOUBLE_LENGTH];

	for (int otherColumnIndex = otherColumnIndexStart; otherColumnIndex <= otherColumnIndexEnd; otherColumnIndex++) {

		for (i = 0; i < commonDimLength; i++) {
			otherColumnPtr[i] = other.rowPtr[i][otherColumnIndex];
		}

		//for evry row in matrix
		for (i = 0; i < this->rowCount; i++) {
			accVector = _mm256_setzero_pd();
			//multiply and add elems 
			for (j = 0; j < commonDimLength - commonDimLength % SIMD_DOUBLE_LENGTH; j += SIMD_DOUBLE_LENGTH) {
				accVector = _mm256_fmadd_pd(_mm256_loadu_pd(&this->rowPtr[i][j]), _mm256_loadu_pd(&otherColumnPtr[j]), accVector);
			}

			_mm256_storeu_pd(accArray, accVector);
			for (k = 1; k < SIMD_DOUBLE_LENGTH; k++) {
				accArray[0] += accArray[k];
			}

			for (j = commonDimLength - commonDimLength % SIMD_DOUBLE_LENGTH; j < commonDimLength; j++) {
				accArray[0] += this->rowPtr[i][j] * otherColumnPtr[j];
			}
			retMatrix.rowPtr[i][otherColumnIndex] = accArray[0];
		}
	}
	delete[] otherColumnPtr;
}

template <>
MyAlgebra::CMatrix<int> MyAlgebra::CMatrix<int>::addMatrixOperation(const CMatrix<int>& other)const {
	if (!equalDimensions(other)) {
		throw DimensionMismatchException(OP_ADD, getDims(), other.getDims());
	}

	CMatrix<int> retMatrix(*this);

	if (!isInitialized(retMatrix)) {
		throw MatrixNotInitialized(OP_ADD);
	}

	int j;
	for (int i = 0; i < rowCount; i++) {
		for (j = 0; j < columnCount - columnCount % SIMD_INT_LENGTH; j += SIMD_INT_LENGTH) {
			_mm256_storeu_epi32(&retMatrix.rowPtr[i][j], _mm256_add_epi32(
				_mm256_loadu_epi32(&retMatrix.rowPtr[i][j]), _mm256_loadu_epi32(&other.rowPtr[i][j])));
		}
		for (; j < columnCount; j++) {
			retMatrix.rowPtr[i][j] += other.rowPtr[i][j];
		}
	}
	return std::move(retMatrix);
}


template <>
MyAlgebra::CMatrix<float> MyAlgebra::CMatrix<float>::addMatrixOperation(const CMatrix<float>& other)const {
	if (!equalDimensions(other)) {
		throw DimensionMismatchException(OP_ADD, getDims(), other.getDims());
	}

	CMatrix<float> retMatrix(*this);

	if (!isInitialized(retMatrix)) {
		throw MatrixNotInitialized(OP_ADD);
	}

	int j;
	for (int i = 0; i < rowCount; i++) {
		for (j = 0; j < columnCount - columnCount % SIMD_FLOAT_LENGTH; j += SIMD_FLOAT_LENGTH) {
			_mm256_storeu_ps(&retMatrix.rowPtr[i][j], _mm256_add_ps(
				_mm256_loadu_ps(&retMatrix.rowPtr[i][j]), _mm256_loadu_ps(&other.rowPtr[i][j])));
		}
		for (; j < columnCount; j++) {
			retMatrix.rowPtr[i][j] += other.rowPtr[i][j];
		}
	}

	return std::move(retMatrix);
}

template <>
MyAlgebra::CMatrix<double> MyAlgebra::CMatrix<double>::addMatrixOperation(const CMatrix<double>& other)const {
	if (!equalDimensions(other)) {
		throw DimensionMismatchException(OP_ADD, getDims(), other.getDims());
	}

	CMatrix<double> retMatrix(*this);

	if (!isInitialized(retMatrix)) {
		throw MatrixNotInitialized(OP_ADD);
	}

	int j;
	for (int i = 0; i < rowCount; i++) {
		for (j = 0; j < columnCount - columnCount % SIMD_DOUBLE_LENGTH; j += SIMD_DOUBLE_LENGTH) {
			_mm256_storeu_pd(&retMatrix.rowPtr[i][j], _mm256_add_pd(
				_mm256_loadu_pd(&retMatrix.rowPtr[i][j]), _mm256_loadu_pd(&other.rowPtr[i][j])));
		}
		for (; j < columnCount; j++) {
			retMatrix.rowPtr[i][j] += other.rowPtr[i][j];
		}
	}
	return std::move(retMatrix);
}

template <>
bool MyAlgebra::CMatrix<int>::comparisionOperation(const  MyAlgebra::CMatrix<int>& other) const {

	if (!equalDimensions(other)) {
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

	if (!equalDimensions(other)) {
		return false;
	}

	for (int i = 0, j = 0; i < this->rowCount; i++) {
		for (j = 0; j < this->columnCount; j++) {
			if (abs(this->rowPtr[i][j] - other.rowPtr[i][j]) > ALG_PRECISION) {
				return false;
			}
		}
	}
	return true;
}

template <>
bool MyAlgebra::CMatrix<double>::comparisionOperation(const  MyAlgebra::CMatrix<double>& other) const {

	if (!equalDimensions(other)) {
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
