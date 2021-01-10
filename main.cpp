
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

#include <stdlib.h>
#include "CMatrix.h"

int main() {

	//MyAlgebra::CMatrix<int> * mat = new MyAlgebra::CMatrix<int>(2,2,true);
	MyAlgebra::CMatrix<float>* mat = new MyAlgebra::CMatrix<float>(2, (float)5.0);

	MyAlgebra::CMatrix<float>* mat2 = new MyAlgebra::CMatrix<float>(*mat);

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
	MyAlgebra::CMatrix<float> kol = mat->getColumnVector(0);
	kol.display();
	MyAlgebra::CMatrix<float> row = mat->getRowVector(1);
	row.display();

	MyAlgebra::CMatrix<float>(3, 3, true).display();

	(mat->getColumnVector(0)).display(); //ask about it

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
		(mat2->multiply(*mat3)).display();
	}

	
	std::cout << TILDES << std::endl;
	{
		mat2->display();
		mat3->display();
		(*mat2 + (*mat3)).display();
		(mat2->add(*mat3)).display();
	}
	
	
	std::cout << TILDES << std::endl;
	{
		mat2->display();
		mat3->display();
		(*mat2 - (*mat3)).display();
		(mat2->substract(*mat3)).display();
	}

	std::cout << TILDES << std::endl;
	{
		mat2->display();
		(-(*mat2)).display();
		(mat2->unary()).display();

	}


	std::cout << TILDES << std::endl;
	delete mat;
	{
		mat = new MyAlgebra::CMatrix<float>(2, 3,true);

		mat->display();
		(~(*mat)).display();
		(mat->transponse()).display();
	}
	
	delete mat;
	delete mat2;
	delete mat3;

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2,true);
		matrix.display();
		matrix.power(0).display();
		matrix.power(1).display();
		matrix.power(2).display();
	}

	std::cout << TILDES << std::endl;
	{
		MyAlgebra::CMatrix<int> matrix = MyAlgebra::CMatrix<int>(2, 2, true);
		MyAlgebra::CMatrix<int> matrix2 = MyAlgebra::CMatrix<int>(2, 2, true);

		matrix.display();
		matrix2.display();

		matrix.dotProduct(matrix2).display();
	}


	_CrtDumpMemoryLeaks();
	

	return 0;
}
