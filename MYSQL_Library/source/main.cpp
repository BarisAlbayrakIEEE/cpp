// main.cpp

#include <iostream>
#include <string>
#include "headerMYSQL.hxx"

using namespace NameMYSQL;

int main()
{
    // Initialize the DB parameters
    const char* hostName = "localhost";
    unsigned int port = 3306;
    const char* DBName = "sakila";
    const char* userName = "local_aee";
    const char* password = "KargaCrow01";

    // Create the MYSQL connection
    MYSQL* connection;
    try {
        connection = NameMYSQL::MYSQL_Library::createConnection(
            hostName,
            port,
            DBName,
            userName,
            password);
    }
    catch (...) {
        return -1;
    }

    // Perform the query
    const char* query = "SELECT last_update FROM actor WHERE first_name = 'PENELOPE' AND last_name = 'GUINESS'";
    MYSQL_RES* resultSet;
    try {
        resultSet = NameMYSQL::MYSQL_Library::createQuery(connection, query);
    }
    catch (...) {
        return -1;
    }

    // Get the data
    MYSQL_ROW row;
    try {
        auto count = (int)mysql_num_rows(resultSet);
        printf("row count: %d\n", count);
        while (true)
        {
            row = mysql_fetch_row(resultSet);
            if (!row) { break; }
            for (int i = 0; i < mysql_num_fields(resultSet); i++)
            {
                printf("%s \t", row[i] != nullptr ? row[i] : "NULL");
            }
            printf("\n");
        }
    }
    catch (...) {
        return -1;
    }

    return 0;
}
