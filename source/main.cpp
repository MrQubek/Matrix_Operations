
#include "simpleTests.h"

#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

#include <iostream>

#include "CMatrix.h"

// W poniższym pliku - impelementacja referencyjna
#include "RefAlgebra/CMtx.h"

#include <stdlib.h>
#include <stdio.h>


// ===================================================================
// FUNKCJE DO POMIARU CZASU
// ===================================================================

#ifdef _WIN32
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif
#include <time.h>
#include <math.h>

double mygettime(void) {
# ifdef _WIN32
	struct _timeb tb;
	_ftime_s(&tb);	// W systemie 64-bitowym może zajść potrzeba zmiany nazwy funkcji na "_ftime64_s"
	return (double)tb.time + (0.001 * (double)tb.millitm);
# else
	struct timeval tv;
	if (gettimeofday(&tv, 0) < 0) {
		perror("oops");
	}
	return (double)tv.tv_sec + (0.000001 * (double)tv.tv_usec);
# endif
}

// ===================================================================
// FUNKCJA OCENY CZASU WYKONANIA
// ===================================================================

// Definiujemy szablon aby łatwiej uruchamiać testy dla róznych implementacji
// klasy. Różne implementacje będą umieszczone w różnych przestrzeniach nazw.
template<typename T>
double test()
{
	// Przykładowe testowe obliczenie macierzowe. Podobne obliczenia będą 
	// używane do oceny efektywności implementacji w konkursie.
	const int SIZE = 1000;
	const int ITER_CNT = 10;

	T A(SIZE, SIZE, true);
	T B(SIZE, SIZE, true);
	T W(1, 1, false);
	double t1 = mygettime();

	for (int i = 0; i < ITER_CNT; i++)
	{
		B = ((0.1f * i) * A + B * B) * 1.e-4f;
		B = -B * ~(A + B);
	}
	W = (B + A);

	double exec_time = mygettime() - t1;

	//W.display();

	return exec_time;
}


int main(int argc, char* argv[])
{
	if(1){
	double t_prog = test<MyAlgebra::CMatrix<float>>();
	double t_ref = 93;//test<RefAlgebra::CMtx>();

	printf("Czas wykonania referencyjny: %7.2lfs\n", t_ref);
	printf("Czas wykonania testowany:    %7.2lfs\n", t_prog);
	printf("Wspolczynnik przyspieszenia Q: %5.2lf", t_ref / t_prog);
	}
	else {
		try{
		simpleTests();
		}
		catch (std::exception & e){
			std::cout << e.what() << std::endl;
		}
	}
	_CrtDumpMemoryLeaks();

	return 0;
}
