#pragma once

#define RANDOM_MAX_VALUE 4
#define RANDOM_MIN_VALUE 0

#define SEPARATOR " "

#define TILDES "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#define EX_NOT_SQUARE_MATRIX "Not square matrix" 
#define EX_MATRIX_NOT_INITIALIZED "rowPtr equals NULL"

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

enum class ErrorCodes : int {
	OK = 0,
	FILE_NOT_OPEN,
	BAD_COLUMN_COUNT,
	BAD_ROW_COUNT,
	BAD_VALUE

};