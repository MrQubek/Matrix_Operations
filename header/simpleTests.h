#pragma once

#include <stdlib.h>
#include "CMatrix.h"

inline void simpleTests() {

	//MyAlgebra::CMatrix<int> * mat = new MyAlgebra::CMatrix<int>(2,2,true);
	MyAlgebra::CMatrix<float>* mat = new MyAlgebra::CMatrix<float>(2, (float)5.0);

	MyAlgebra::CMatrix<float>* mat2 = new MyAlgebra::CMatrix<float>(*mat);

	MatrixError errCode = MatrixError();

	mat2->set(0, 0, 4.0);

	mat->display();
	mat2->display();

	std::cout << TILDES << std::endl;

	MyAlgebra::CMatrix<float>* mat3 = new MyAlgebra::CMatrix<float>(std::move(*mat2));

	mat->display();
	mat2->display();
	std::cout << std::endl;
	mat3->display();

	std::cout << TILDES << std::endl;

	*mat2 = std::move(*mat);

	mat->display();
	std::cout << std::endl;
	mat2->display();
	mat3->display();

	std::cout << TILDES << std::endl;

	*mat3 = *mat2;
	mat->display();
	mat2->display();
	mat3->display();


	std::cout << TILDES << std::endl;

	mat3->toDiagonal(4.5).display();
	(*mat3 = 1.5).display();
	std::cout << TILDES << std::endl;


	{
		MyAlgebra::CMatrix<int> kol(2, 2, true);
		kol.getColumnVector(0, errCode).display();

		kol.getColumnVector(2, errCode).display();
		if(errCode.getErrorCode() == MatrixErrorCode::WRONG_ARGUMENT) {
			std::cout << errCode.getOpName() << NEW_LINE;
		}

		MyAlgebra::CMatrix<int> row(2, 2, true);
		row.getRowVector(1, errCode).display();

	}

	std::cout << TILDES << std::endl;
	{
		(*mat3 * 4.5).display();
	}

	std::cout << TILDES << std::endl;
	{
		mat2->display();
		mat3->display();
		MyAlgebra::CMatrix<float> matMultiply = (*mat2 * (*mat3));
		matMultiply.display();
		(mat2->multiply(*mat3, errCode)).display();
	}


	std::cout << TILDES << std::endl;
	{
		mat2->display();
		mat3->display();
		(*mat2 + (*mat3)).display();
		(mat2->add(*mat3, errCode)).display();
	}


	std::cout << TILDES << std::endl;
	{
		mat2->display();
		mat3->display();
		(*mat2 - (*mat3)).display();
		(mat2->substract(*mat3, errCode)).display();
	}

	std::cout << TILDES << std::endl;
	{
		mat2->display();
		(-(*mat2)).display();
		(mat2->unary(errCode)).display();

	}

	std::cout << TILDES << std::endl;
	delete mat;
	{
		mat = new MyAlgebra::CMatrix<float>(2, 3, true);

		mat->display();
		(~(*mat)).display();
		(mat->transponse(errCode)).display();
	}

	delete mat;
	delete mat2;
	delete mat3;

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, true);
		matrix.display();
		matrix.power(0, errCode).display();
		matrix.power(1, errCode).display();
		matrix.power(2, errCode).display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, true);
		MyAlgebra::CMatrix<int> matrix2 = MyAlgebra::CMatrix<int>(2, 2, true);

		matrix.display();
		matrix2.display();

		matrix.dotProduct(matrix2, errCode).display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<float> matrix = MyAlgebra::CMatrix<float>(2, 2, true);
		MyAlgebra::CMatrix<float> matrix2 = MyAlgebra::CMatrix<float>(2, 2, true);
		MyAlgebra::CMatrix<float> matrix3 = MyAlgebra::CMatrix<float>(3, 2, true);

		std::cout << (0.0f < 10e-6f) << std::endl;
		std::cout << (matrix == matrix) << std::endl;
	}

	
	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<float> matrix = MyAlgebra::CMatrix<float>(2, 2, false);

		if(matrix.readMatrixFromFile("floatMatrix.txt"))
		matrix.display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, false);

		if (matrix.readMatrixFromFile("intMatrix.txt"))
			matrix.display();
	}
	
	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, false);

		if (matrix.readMatrixFromFile("badIntMatrix_1.txt"))
			matrix.display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, false);

		if (matrix.readMatrixFromFile("badIntMatrix_1.txt"))
			matrix.display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, false);

		matrix[0][0] = 0;
		matrix[0][1] = 1;
		matrix[1][0] = 2;
		matrix[1][1] = 3;

		matrix.display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, true);
		MyAlgebra::CMatrix<int> matrix2 = MyAlgebra::CMatrix<int>(2, 2, true);
		MyAlgebra::CMatrix<int> matrix3 = MyAlgebra::CMatrix<int>(2, 2, true);

		matrix.display();
		matrix2.display();
		matrix3.display();

		(matrix* matrix2).display();

		matrix = matrix * matrix2 - matrix3;

		matrix.display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, true);
		MyAlgebra::CMatrix<int> matrix2 = MyAlgebra::CMatrix<int>(2, 2, true);
		MyAlgebra::CMatrix<int> matrix3 = MyAlgebra::CMatrix<int>(2, 2, true);

		matrix.display();
		matrix2.display();
		matrix3.display();

		matrix = matrix.dotProduct(matrix2,errCode) - matrix3;

		if (errCode) {
			std::cout << OKK << NEW_LINE;
		} else {
			std::cout << ERROR << NEW_LINE;
		}

		matrix.display();
	}

}