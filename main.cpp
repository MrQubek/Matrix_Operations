
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
	(*mat3 = 2.5).display();
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

	delete mat;
	delete mat2;
	delete mat3;

	_CrtDumpMemoryLeaks();


	return 0;
}