
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

	*mat = *mat2;
	mat->display();
	mat2->display();
	mat3->display();

	delete mat;
	delete mat2;
	delete mat3;

	_CrtDumpMemoryLeaks();


	return 0;
}