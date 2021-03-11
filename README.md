# MatrixOperations
Template class which allows creation of matrixes and operations on them. 

### Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Details about implementation](#details-about-implementation)
* [Description of header files](#description-of-header-files)
* [Compilation](#compilation)
* [Sources](#sources)


### General info
The implementation uses threads and Intel AVX/AVX2 SIMD instructions in order to speed up calculation time. Move semantics is widely used in the project as well.

Initially, the project was created as a part of passing the course of Techniques of Effective Programming (pol. Techinik Efektywnego Programowania). Multithreading and SIMD instructions have been added later in order to learn basics how to use them.

### Technologies
The project is created with:
* Visual Studio 2019
* Intel AVX/AVX2 instructions
* Move semantics
* C++11 Threads
* Templates

### Details about implementation
Matrixes and operations (through class _CMatrix_) are defined for C++ types int, float and double. 

The constructor _CMatrix(int rowCnt, int colCnt, bool randInit)_ initialises matrix by allocating a memory for it. If the "randInit" variable is true, the cells of matrix will be filled with uniformly distributed random numbers constrained by constants defined in header const.h; otherwise the cells will be left with unspecified values.

The followig operations may be performed on the matrixes:
* multiplication of matrixes
* multiplication by constant
* addition
* substarction
* unary operation
* transposition
* exponentiation
* dot product of two matrixes (defined as A<sup>T</sup>*B)

Most operations are available through either functions or overloaded operators. When calling an operation by a function, the user has to pass a MatrixError object to a function. In case of an error (e.g. dimensions mismatch of matrixes), a function will return a zero matrix and it will put an error code into the MatrixError object (both defined in header MatrixError.h), whereas an overloaded operator will throw a custom exception (defined in header MatrixException.h).

The call of an operation causes new memory allocation for its result and returns the result via move sematic. 

Most of the operations are based on other operations implemented with Intel AVX/AVX2 SIMD instructions. Those instructions are used in functions _thread\_OpName\_Operation_, where _\_OpName\__ is one of the following:
* MultiplyConstant
* Multiply
* Add
* Substract

Using an operation usually causes the following call tree: OpName (optional) -> overloaded operator -> OpNameOperation -> threadOpNameOperation.

### Content of header files
##### CMatrix.h
CMatrix temaplate class is declared here. 

##### const.h
Various constant values used in the project are defined here. A user can use them to specify among others:
* maximum thread count;
* threshold when to use multiple threads when adding, substracting or multipling matrix by constant value;
* maximum differential of float/double values when comparing two matrixes;
* range of values when filling matrixes with random numbers;
* length of SIMD vectors.  

##### FileError.h
FileError class and error codes dedicated malfunctions during reading matrix from file are defined here.

##### MatrixError.h
MatrixError class and error codes dedicated malfunctions during mathematical operations on matrixes are defined here. The error code is returned when an operation is called by a function. 


In case of no malfunction present during an operation, the returned error code would be 

##### MatrixException.h
Classes inheriting from std::exception dedicated malfunctions during mathematical operations on matrixes are defined here. An exception is thrown when an operation is called by overloaded operator.

### Compilation
A benchmark has been prepeared in order to show the speed of calculation time of operations.  Build the project with MSBuild by executing command below in cmd:

```
MSBuild /nologo /v:n /p:Configuration=Release
```
and then execute compiled Matrix_Operations.exe file.

### Sources
Minimal interface of Class CMatrix has been imposed by Dr Michał Przewoźniczek, which was later expanded and impemented by myself.
