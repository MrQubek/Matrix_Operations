#pragma once

#include <stdlib.h>
#include "CMatrix.h"

#define TYPE int

inline void simpleTests() {

	MatrixError errCode = MatrixError();

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<TYPE> matrix = MyAlgebra::CMatrix<TYPE>(4, 27, true);
		MyAlgebra::CMatrix<TYPE> matrix2 = MyAlgebra::CMatrix<TYPE>(27, 6, true);


		matrix.display();
		(~matrix2).display();
		matrix2.display();


		std::cout << NEW_LINE;

		(matrix + matrix2).display();
		std::cout << NEW_LINE;

		(matrix + matrix2).display();

		if (errCode) {
			std::cout << OKK << NEW_LINE;
		} else {
			std::cout << ERROR << NEW_LINE;
		}
	}

}