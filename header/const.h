#pragma once

// =========================================================================
// control constants
// =========================================================================

#define MAX_THREAD_COUNT 4
#define MIN_CELLS_FOR_MULTITHREADING_IN_SIMPLE_OPS 25

#define MAX_COMPARE_ERROR 1e-6f

// 0 - no limits on value of generated random numbers
// 1 - value of generated random numbers constrained by RANDOM_MIN_VALUE and RANDOM_MAX_VALUE
#define LIMIT_MAX_RANDOM_VALUE 1
#define RANDOM_MAX_VALUE 2
#define RANDOM_MIN_VALUE 0

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