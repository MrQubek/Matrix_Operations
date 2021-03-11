

#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

#include <chrono>

#include "CMatrix.h"
#include "simpleTests.h"

// template function for showing speed of matrix operations.
// description of constats:
//		SIZE - size of matrixes.
//		ITER_CNT - number of repetitions of operations in order get average operation time.
//		DISPLAY_RESULTS - display matrix resulted from matrix operation.
//		CONSTANT - consatnt used in multiplication by constant operation.
template<typename T>
void benchmark() {

	const int SIZE = 2001;
	const int ITER_CNT = 5;
	const bool DISPLAY_RESULTS = false;
	const int CONSTANT = 2;

	T A(SIZE, SIZE, true);
	T B(SIZE, SIZE, true);
	T C(1, 1, false);

	try {
		std::cout << "Size of matrixes: " << SIZE << "x" << SIZE << std::endl;

		if (DISPLAY_RESULTS) {
			std::cout << "Matrix A:\n";
			A.display();
			std::cout << "Matrix B:\n";
			B.display();
		}

		std::cout << "Executing multiplying matrixes operation C = A * B."<<std::endl;
		auto operationBeginTime = std::chrono::steady_clock::now();
		for (int i = 0; i < ITER_CNT;i++) {
			C = A * B;
		}
		auto operationEndTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> totalOperationTime = operationEndTime - operationBeginTime;
		std::cout << "Average multiplication time: "<< totalOperationTime.count()/ITER_CNT << "s.\n";

		if (DISPLAY_RESULTS) {
			std::cout << "\nMatrix C\n";
			C.display();
			std::cout << std::endl;
		}

		std::cout << "Executing adding matrixes operation C = A + B." << std::endl;
		operationBeginTime = std::chrono::steady_clock::now();
		for (int i = 0; i < ITER_CNT; i++) {
			C = A + B;
		}
		operationEndTime = std::chrono::steady_clock::now();
		totalOperationTime = operationEndTime - operationBeginTime;
		std::cout << "Average addition time: " << totalOperationTime.count() / ITER_CNT << "s.\n";

		if (DISPLAY_RESULTS) {
			std::cout << "\nMatrix C\n";
			C.display();
			std::cout << std::endl;
		}

		std::cout << "Executing substracting matrixes operation C = A - B." << std::endl;
		operationBeginTime = std::chrono::steady_clock::now();
		for (int i = 0; i < ITER_CNT; i++) {
			C = A - B;
		}
		operationEndTime = std::chrono::steady_clock::now();
		totalOperationTime = operationEndTime - operationBeginTime;
		std::cout << "Average substraction time: " << totalOperationTime.count() / ITER_CNT << "s.\n";

		if (DISPLAY_RESULTS) {
			std::cout << "\nMatrix C\n";
			C.display();
			std::cout << std::endl;
		}

		std::cout << "Executing multiplication by constant operation C = A * CONSTANT." << std::endl;
		operationBeginTime = std::chrono::steady_clock::now();
		for (int i = 0; i < ITER_CNT; i++) {
			C = A * CONSTANT;
		}
		operationEndTime = std::chrono::steady_clock::now();
		totalOperationTime = operationEndTime - operationBeginTime;
		std::cout << "Average  multiplication by constant time: " << totalOperationTime.count() / ITER_CNT << "s.\n";

		if (DISPLAY_RESULTS) {
			std::cout << "\nMatrix C\n";
			C.display();
			std::cout << std::endl;
		}

		std::cout << "Press any key to exit.\n";
		std::cin.get();
	}
	catch (std::exception& e) {
		std::cout << "Exception captured in: " << e.what() << "\nAborting benchmark\n";
	}
}


int main(int argc, char* argv[])
{

	benchmark< MyAlgebra::CMatrix<int>>();

	_CrtDumpMemoryLeaks();

	return 0;
}
