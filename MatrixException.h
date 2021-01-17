#pragma once

#include <exception>
#include <string>

class MatrixException : public std::exception {
protected:
	std::string operationName;
public:
	MatrixException(std::string opName)
		: operationName(opName) {}

	std::string getOpName() const { return operationName; }
};

class MatrixNotInitialized : public MatrixException {
public:
	MatrixNotInitialized(std::string opName)
		: MatrixException(opName) {}
};

class WrongArgument : public MatrixException {
protected:
	int wrongArgument1;
public:
	WrongArgument(std::string opName, int wrArg1)
		: MatrixException(opName), wrongArgument1(wrArg1) {}

	int getWronArg1() const{ return wrongArgument1; }
};

class TwoWrongArguments: public WrongArgument {
protected:
	int wrongArgument2;
public:
	TwoWrongArguments(std::string opName, int wrArg1, int wrArg2)
		: WrongArgument(opName, wrArg1), wrongArgument2(wrArg2) {}
	int getWronArg2() const { return wrongArgument2; }
};

class NotSquareMatrix :public MatrixException {
protected:
	std::pair<int, int> dimensions;
public:
	NotSquareMatrix(std::string opName, std::pair<int, int> dim)
		: MatrixException(opName), dimensions(dim) {}
	std::pair<int, int> getDims() const { return dimensions; }
};

class DimensionMismatchException :public NotSquareMatrix {
protected:
	std::pair<int, int> dimensionsOther;
public:
	DimensionMismatchException(std::string opName, std::pair<int, int> dims, std::pair<int, int> dimOther)
		: NotSquareMatrix(opName, dims), dimensionsOther(dimOther) {}

	std::pair<int, int> getDimsOther() const { return dimensionsOther; }

};

