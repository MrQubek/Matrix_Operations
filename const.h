#pragma once

#define RANDOM_MAX_VALUE 4
#define RANDOM_MIN_VALUE 0

#define SEPARATOR " "

#define TILDES "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#define OP_DIAGONAL "Diagonal operation"
#define OP_MULTIPLY_CONST "Multiply by constant operation"
#define OP_MULTIPLY_MATRIX "Multiply matrix operation"
#define OP_ADD "Add operation"
#define OP_SUBSTRACT "Substract operation"
#define OP_UNARY "Unary operation"
#define OP_TRANSPONSE "Transponse operation"
#define OP_POWER "Power operation"
#define OP_DOT_PRODUCT "Dot product operation"
#define OP_COMPARE "Compare operation"

#define OP_GET_ROW "Get row vector operation"
#define OP_GET_COLUMN "Get column vector operation"

#define OP_GET "Get operation"
#define OP_SET "Set operation"


enum class FileErrorCode : int {
	OK = 0,
	FILE_NOT_OPEN,
	BAD_FIRST_GETLINE,
	BAD_COLUMN_COUNT,
	BAD_MEM_ALLOC,
	BAD_GETLINE,
	READ_BAD_VALUE
};

enum class MatrixErrorCode : int {
	OK = 0,
	MATRIX_NOT_INITIALIZED,
	WRONG_ARGUMENT,
	TWO_WRONG_ARGUMENTS,
	NOT_SQUARE_MATRIX,
	DIMENSION_MISMATCH
};