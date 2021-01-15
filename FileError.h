#pragma once

#include "const.h"

class FileError
{
private:
	ErrorCodes errorCode;
public:
	FileError(ErrorCodes eCode) :errorCode(eCode) {}

	ErrorCodes getErrorCode() { return errorCode; }

	operator bool() const { return errorCode == ErrorCodes::OK; }
};

