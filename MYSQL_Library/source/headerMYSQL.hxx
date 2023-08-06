/// <summary>
/// This project contains MYSQL methods
/// </summary>

#pragma warning(disable : 4290)

#ifndef _NameMYSQL_HeaderFile
#define _NameMYSQL_HeaderFile

#include <iostream>
#include <time.h>
#include <stdio.h>
#include "mysql.h"

namespace NameMYSQL {
	class MYSQL_Library
	{
	public:
		static MYSQL* createConnection(
			const char* theHost,
			unsigned int thePort,
			const char* theDBName,
			const char* theUserName,
			const char* thePassword);
		static MYSQL_RES* createQuery(
			const char* theHost,
			unsigned int thePort,
			const char* theDBName,
			const char* theUserName,
			const char* thePassword,
			const char* theQuery);
		static MYSQL_RES* createQuery(
			MYSQL* theConnection,
			const char* theQuery);
	};
}

#endif
