#pragma once

#include "const.h"

class FileError
{
private:
	FileErrorCode errorCode;
public:
	FileError(FileErrorCode eCode) :errorCode(eCode) {}

	FileErrorCode getErrorCode() { return errorCode; }

	operator bool() const { return errorCode == FileErrorCode::OK; }
};

