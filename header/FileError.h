#pragma once

#include "const.h"

enum class FileErrorCode : int {
	OK = 0,
	FILE_NOT_OPEN,
	BAD_FIRST_GETLINE,
	BAD_COLUMN_COUNT,
	BAD_MEM_ALLOC,
	BAD_GETLINE,
	READ_BAD_VALUE
};

class FileError
{
private:
	FileErrorCode errorCode;
public:
	FileError(FileErrorCode eCode) : errorCode(eCode) {}

	FileErrorCode getErrorCode() { return errorCode; }

	operator bool() const { return errorCode == FileErrorCode::OK; }
};

