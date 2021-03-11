#pragma once

// =========================================================================
// control constants
// =========================================================================

#define MAX_THREAD_COUNT 4

// number of cells in matrix is defined as rowCount * columnCount
#define MIN_CELLS_FOR_MULTITHREADING_IN_SIMPLE_OPS 25000000

#define MAX_COMPARE_ERROR 1e-6f

// 1 - value of generated random numbers constrained by CUSTOM_RANDOM_MIN_VALUE and CUSTOM_RANDOM_MAX_VALUE.
// 0 - default limits on value of generated random numbers:
//		int - [std::numeric_limits<int>::min(), std::numeric_limits<int>::max()]
//		float - [-FLT_MAX/2, FLT_MAX/2]
//		double - [-DBL_MAX/2, DBL_MAX/2]
#define CUSTOM_LIMITS_ON_RANDOM_VALUES 0

#define CUSTOM_RANDOM_MAX_VALUE 1000
#define CUSTOM_RANDOM_MIN_VALUE -1000

// length of SIMD vectors
#define SIMD_INT_LENGTH 8
#define SIMD_FLOAT_LENGTH 8
#define SIMD_DOUBLE_LENGTH 4

// =========================================================================
// string constants
// =========================================================================

#define SEPARATOR " "
#define NEW_LINE "\n"
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

#define OP_THREAD " thread"

#define OKK "OKK"
#define ERROR "ERROR"