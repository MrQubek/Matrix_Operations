#pragma once

#include<utility>
#include<string>
#include<cstdlib>
#include<time.h>
#include <random>
#include <iostream>

#include "const.h"


namespace MyAlgebra
{
	template <typename T>
	class CMatrix
	{
	private:
		T** rowPtr;

		int     rowCount;
		int     columnCount;

		// return:	true if memory allocation for matrix successful
		//			false otherwise
		// warning: possible memory leaks if memory was already allocated
		bool allocateMemory(int rowCnt, int colCnt);

		void deallocateMemory();

		void populateMatrixWithRandomNumbers();

		void copyMatrixValues(const CMatrix& other);

		void copyOperation(const CMatrix& other);

		void moveOperation(CMatrix&& other);

		void diagonalOperation(T diagonal);

		CMatrix<T>& multiplyConstantOperation(T multiplier);

		CMatrix<T> multiplyMatrixOperation(const CMatrix<T> & other);

		CMatrix(T** newRowPtr, int rowCnt, int colCnt);

	public:
		static const float ALG_PRECISION;

		// =========================================================================
		// Constructors:
		// =========================================================================

		// create matrix with (possible) random values
		CMatrix(int rowCnt, int colCnt, bool randInit = false);

		// create square diagonal matrix
		CMatrix(int rowCnt, T diagonal);

		CMatrix(const CMatrix& other);

		CMatrix(CMatrix&& other) noexcept;

		~CMatrix();

		CMatrix createMatrix(int rowCnt, int colCnt, bool randInit = false) {
			return std::move(CMatrix(rowCnt, colCnt, randInit));
		}

		// =========================================================================
		// Assign operators:
		// =========================================================================

		const CMatrix& operator=(const CMatrix& other);
		const CMatrix& assignMatrix(const CMatrix& other);

		// move operator
		const CMatrix& operator=(CMatrix&& other)noexcept;
		const CMatrix& moveMatrix(CMatrix&& other);

		// switch matrix to diagonal matrix 
		const CMatrix& operator=(T diagonal);
		const CMatrix& toDiagonal(T diagonal);

		// =========================================================================
		// indexing matrix
		// =========================================================================

		T* operator[](int row_ind) { return rowPtr[row_ind]; }
		T* getRow(int row_ind) { return rowPtr[row_ind]; }

		void set(int x, int y, T newVal) {
			if (rowPtr != nullptr)rowPtr[x][y] = newVal;
		}

		T get(int x, int y) {
			//if (rowPtr != nullptr && x < columnCount && y < rowCount)
			return rowPtr[x][y];
		}

		int getRowCount() { return rowCount; }
		int getColumnCount() { return columnCount; }

		// =========================================================================
		// return vector matrix from rows/columns
		// =========================================================================

		CMatrix getRowVector(int rowIndex);

		CMatrix getColumnVector(int columnIndex);

		// =========================================================================
		// algebraic operations
		// =========================================================================

		CMatrix operator*(const CMatrix& other) const;
		CMatrix multiply(const CMatrix& other) const;

		// multiply matrix by constant
		CMatrix operator*(T multiplier) const;
		CMatrix multiply(T multiplier) const;


		CMatrix operator+(const CMatrix& other) const;
		CMatrix add(const CMatrix& other) const;

		CMatrix operator-(const CMatrix& other) const;
		CMatrix substract(const CMatrix& other) const;


		// change sign of all elements of matrix
		CMatrix operator-() const;
		CMatrix unary() const;

		// transponse matrix
		CMatrix operator~() const;
		CMatrix transponse() const;

		// accept only power >= 0:
		//    power = 0  - return unity matrix
		//    power = 1  - return copy of matrix
		//    power > 1  - return power of matrix
		CMatrix operator^(int power) const;
		CMatrix power(int power) const;

		// dot product (iloczyn skalarny) A^T*B
		CMatrix dotProduct(const CMatrix& other) const;

		// compare matrix with accuracy to ALG_PRECISION
		bool operator==(const CMatrix& rhs) const;
		bool compareTo(const CMatrix& rhs) const;

		// =========================================================================
		// algebraic operations
		// =========================================================================

		bool readMatrixFromFile(std::string fileName);

		// only for tests - display matrix by rows on stdout
		void display() const;

		// friend CMtx operator*( T multiplier, const CMtx &rhs );
	};

	//template <typename T>
	//CMatrix<T> operator*(T multiplier, const CMatrix<T>& rhs);

	// =========================================================================
	// methods declarations
	// =========================================================================

	template <typename T>
	bool CMatrix<T>::allocateMemory(int rowCnt, int colCnt) {
		//to do: add exceprions
		if (rowCnt < 0 || colCnt < 0) {
			rowPtr = nullptr;
			return false;
		}

		rowPtr = new(std::nothrow) T * [rowCnt];

		if (rowPtr != nullptr) {
			for (int i = 0; i < rowCnt; i++) {
				rowPtr[i] = new(std::nothrow) T[colCnt];

				if (rowPtr[i] == nullptr) {
					for (; i > 0; i--) {
						delete[] rowPtr[i];
					}
					delete[] rowPtr;
					rowPtr = nullptr;
					return false;
				}
			}
		}
		return true;
	}

	template <typename T>
	void CMatrix<T>::deallocateMemory() {
		for (int i = 0; i < rowCount; i++) {
			delete[] rowPtr[i];
		}
		delete[] rowPtr;
	}

	template <typename T>
	void CMatrix<T>::copyMatrixValues(const CMatrix<T>& other) {
		for (int i = 0, j = 0; i < other.columnCount; i++) {
			for (j = 0; j < other.columnCount; j++) {
				this->rowPtr[i][j] = other.rowPtr[i][j];
			}
		}
	}

	template <typename T>
	void CMatrix<T>::copyOperation(const CMatrix& other) {
		if (allocateMemory(other.rowCount, other.columnCount)) {
			this->rowCount = other.rowCount;
			this->columnCount = other.columnCount;

			copyMatrixValues(other);

		}
		else {
			this->rowCount = 0;
			this->columnCount = 0;
		}
	}

	template <typename T>
	void CMatrix<T>::moveOperation(CMatrix&& other) {

		this->rowPtr = other.rowPtr;
		this->rowCount = other.rowCount;
		this->columnCount = other.columnCount;

		other.rowPtr = nullptr;
		other.rowCount = 0;
		other.columnCount = 0;

	}

	template <typename T>
	void CMatrix<T>::diagonalOperation(T diagonal) {

		for (int i = 0, j = 0; i < rowCount; i++) {
			for (j = 0; j < columnCount; j++) {
				if (i == j) {
					rowPtr[i][j] = diagonal;
				}
				else {
					rowPtr[i][j] = 0;
				}
			}
		}
	}
	template <typename T>
	CMatrix<T>& CMatrix<T>::multiplyConstantOperation(T multiplier) {
		for (int i = 0, j; i < rowCount; i++) {
			for (j = 0; j < columnCount; j++) {
				rowPtr[i][j] *= multiplier;
			}
		}
		return *this;
	}

	/*template <typename T>
	CMatrix<T>& CMatrix<T>::multiplyMatrixOperation(const CMatrix<T>& other) {

		//to do: add exceptions

		CMatrix<T> retMatrix(this->rowCount, other.columnCount);

		if (retMatrix.getRowCount() > 0 && retMatrix.getColumnCount() > 0) {
			for (int i = 0, j, r, fieldSum; i < retMatrix.getRowCount; i++) {
				for (j = 0; j < retMatrix.getColumnCount; j++) {
					for (r = 0, fieldSum = 0; r < this->columnCount; r++) {
						fieldSum += rowPtr[i][r] * other.rowPtr[r][j];
					}
					retMatrix[i][j] = fieldSum;
				}
			}
		}
		return std::move(retMatrix);
	}*/
	
	template <typename T>
	CMatrix<T>::CMatrix(T** newRowPtr, int rowCnt, int colCnt) {
		//to do: add exceptions

		rowPtr = newRowPtr;
		rowCount = rowCnt;
		columnCount = colCnt;
	}

	template <typename T>
	CMatrix<T>::CMatrix(int rowCnt, int colCnt, bool randInit) {
		//to do: add exceptions

		if (allocateMemory(rowCnt, colCnt)) {
			this->rowCount = rowCnt;
			this->columnCount = colCnt;

			if (randInit) {
				srand((unsigned int)time(0));
				populateMatrixWithRandomNumbers();
			}
		}
		else {
			this->rowCount = 0;
			this->columnCount = 0;
		}
	}

	template <typename T>
	CMatrix<T>::CMatrix(int rowCnt, T diagonal) {
		//to do: add exceptions

		if (allocateMemory(rowCnt, rowCnt)) {
			this->rowCount = rowCnt;
			this->columnCount = rowCnt;

			for (int i = 0, j = 0; i < rowCnt; i++) {
				for (j = 0; j < rowCnt; j++) {
					if (i == j) {
						rowPtr[i][j] = diagonal;
					}
					else {
						rowPtr[i][j] = 0;
					}
				}

			}
		}
		else {
			this->rowCount = 0;
			this->columnCount = 0;
		}
	}

	template <typename T>
	CMatrix<T>::CMatrix(const CMatrix& other) {
		copyOperation(other);
	}

	template <typename T>
	CMatrix<T>::CMatrix(CMatrix&& other) noexcept {
		moveOperation(std::forward<CMatrix>(other));
	}

	template <typename T>
	CMatrix<T>::~CMatrix() {
		deallocateMemory();
	}

	template <typename T>
	const CMatrix<T>& CMatrix<T>::operator=(const CMatrix& other) {
		if (this != &other) {
			if (this->rowPtr != nullptr) {
				deallocateMemory();
			}
			copyOperation(other);
		}
		return *this;
	}

	template <typename T>
	const CMatrix<T>& CMatrix<T>::assignMatrix(const CMatrix& other) {
		if (this != &other) {
			if (this->rowPtr != nullptr) {
				deallocateMemory();
			}
			copyOperation(other);
		}
		return *this;
	}

	template <typename T>
	const CMatrix<T>& CMatrix<T>::operator=(CMatrix<T>&& other) noexcept {

		if (this != &other) {
			if (this->rowPtr != nullptr) {
				deallocateMemory();
			}
			moveOperation(std::forward<CMatrix>(other));
		}
		return *this;
	}

	template <typename T>
	const CMatrix<T>& CMatrix<T>::moveMatrix(CMatrix<T>&& other) {
		if (this != &other) {
			if (this->rowPtr != nullptr) {
				deallocateMemory();
			}
			moveOperation(std::forward<CMatrix>(other));
		}
		return *this;
	}

	template <typename T>
	const CMatrix<T>& CMatrix<T>::operator=(T diagonal) {
		diagonalOperation(diagonal);
		return *this;
	}

	template <typename T>
	const CMatrix<T>& CMatrix<T>::toDiagonal(T diagonal) {
		diagonalOperation(diagonal);
		return *this;
	}

	template <typename T>
	CMatrix<T> CMatrix<T>::getRowVector(int rowIndex) {
		//to do: add exceptions

		T** newRowPtr = new(std::nothrow) T * [1];
		if (newRowPtr != nullptr) {
			newRowPtr[0] = new(std::nothrow) T[columnCount];

			if (newRowPtr[0] != nullptr) {
				CMatrix retMatrix(newRowPtr, 1, columnCount);
				for (int i = 0; i < columnCount; i++) {
					retMatrix.set(0, i, this->rowPtr[rowIndex][i]);
				}
				return std::move(retMatrix);
			}
			else {
				delete[] newRowPtr;
			}
		}
		return std::move(CMatrix(nullptr, 0, 0));
	}

	template <typename T>
	CMatrix<T> CMatrix<T>::getColumnVector(int columnIndex) {
		//to do: add exceptions

		T** newColumnPtr = new(std::nothrow) T * [rowCount];
		if (newColumnPtr != nullptr) {
			int i = 0;
			for (; i < rowCount; i++)
			{
				newColumnPtr[i] = new T[1];
				if (newColumnPtr[i] == nullptr) {
					for (; i > 0; i--) {
						delete[] newColumnPtr[i];
					}
					delete[] newColumnPtr;
					return std::move(CMatrix(nullptr, 0, 0));
				}
			}

			CMatrix retMatrix(newColumnPtr, rowCount, 1);
			for (i = 0; i < rowCount; i++)
			{
				retMatrix.set(i, 0, this->rowPtr[i][columnIndex]);
			}
			return std::move(retMatrix);
		}
		else {
			return std::move(CMatrix(nullptr, 0, 0));
		}
	}

	template <typename T>
	CMatrix<T> CMatrix<T>::operator*(const CMatrix& other) const {
		//to do: add exceptions

		CMatrix<T> retMatrix(this->rowCount, other.columnCount);

		if (retMatrix.getRowCount() > 0 && retMatrix.getColumnCount() > 0) {
			T fieldSum;
			for (int i = 0, j, r; i < retMatrix.getRowCount(); i++) {
				for (j = 0; j < retMatrix.getColumnCount(); j++) {
					for (r = 0, fieldSum = 0; r < this->columnCount; r++) {
						fieldSum += + rowPtr[i][r] * other.rowPtr[r][j];
					}
					retMatrix[i][j] = fieldSum;
				}
			}
		}
		return std::move(retMatrix);
	}
	template <typename T>
	CMatrix<T> CMatrix<T>::multiply(const CMatrix& other) const {
		return std::move(*this*other);
	}

	template <typename T>
	CMatrix<T> CMatrix<T>::operator*(T multiplier) const {
		return std::move(CMatrix(*this).multiplyConstantOperation(multiplier));
	}
	template <typename T>
	CMatrix<T> CMatrix<T>::multiply(T multiplier) const {
		return std::move(CMatrix(*this).multiplyConstantOperation(multiplier));
	}



	template <typename T>
	void CMatrix<T>::display() const {
		for (int i = 0, j = 0; i < rowCount; i++) {
			for (j = 0; j < columnCount; j++) {
				std::cout << rowPtr[i][j] << SEPARATOR;
			}
			std::cout << std::endl;
		}
	}











}



