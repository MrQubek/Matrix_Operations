#pragma once
#include <string>

enum class MatrixErrorCode : int {
	OK = 0,
	MATRIX_NOT_INITIALIZED,
	WRONG_ARGUMENT,
	TWO_WRONG_ARGUMENTS,
	NOT_SQUARE_MATRIX,
	DIMENSION_MISMATCH,
	BAD_MEM_ALLOC,
	UNKNOWN_ERROR
};

class MatrixError {
private:
	std::string operationName;
	MatrixErrorCode errorCode;

public:
	MatrixError() {
		operationName = "";
		errorCode = MatrixErrorCode::OK;
	}

	MatrixError( std::string opName, MatrixErrorCode eCode) {
		operationName = opName;
		errorCode = eCode;
	}

	std::string getOpName() const {return operationName;}

	MatrixErrorCode getErrorCode() const { return errorCode; }

	operator bool() const {
		return errorCode == MatrixErrorCode::OK;
	}


};