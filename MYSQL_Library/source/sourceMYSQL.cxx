// baris.albayrak.ieee@gmail.com

#include "headerMYSQL.hxx"

namespace NameMYSQL {
    /// <summary>
    /// Creates a MYSQL connection
    /// </summary>
    MYSQL* MYSQL_Library::createConnection(
        const char* theHost,
        unsigned int thePort,
        const char* theDBName,
        const char* theUserName,
        const char* thePassword)
    {
        // Initialize the MYSQL connection
        MYSQL* outConnection;
        outConnection = mysql_init(nullptr);

        // Create the MYSQL connection
        outConnection = mysql_real_connect(
            outConnection,
            theHost,
            theUserName,
            thePassword,
            theDBName,
            thePort,
            nullptr,
            0);
        
        // Inspect whether the connection is successful
        if (!outConnection)
        {
            std::cout << "Failed to connect to the database.";
            return nullptr;
        }

        std::cout << "Successfully connected to the database.";
        return outConnection;
    }

    /// <summary>
    /// Creates a MYSQL connection and performs a MYSQL query
    /// CAUTION:
    ///    Use this method if the connection is intended to be used only once.
    ///    as this method
    ///        creates a connection,
    ///        executes the query and
    ///        deletes the connection.
    ///    Hence, an SQL connection is ccreated and deleted at the end
    ///    everytime this method is called.
    ///    Create a connection and execute queries
    ///    if a number of queries will be executed
    ///    in order for the runtime.
    /// </summary>
    MYSQL_RES* MYSQL_Library::createQuery(
        const char* theHost,
        unsigned int thePort,
        const char* theDBName,
        const char* theUserName,
        const char* thePassword,
        const char* theQuery)
    {
        // Create the MYSQL connection
        MYSQL* connection = MYSQL_Library::createConnection(
            theHost,
            thePort,
            theDBName,
            theUserName,
            thePassword);
        if (!connection) { return nullptr; }

        // Perform the query
        int status = mysql_query(connection, theQuery);
        if (status) {
            std::cout << "Failed to perform the query.";
            return nullptr;
        }
        return mysql_store_result(connection);

        // Close the connection
        mysql_close(connection);
    }

    /// <summary>
    /// Performs a MYSQL query
    /// </summary>
    MYSQL_RES* MYSQL_Library::createQuery(
        MYSQL* theConnection,
        const char* theQuery)
    {
        // Perform the query
        int status = mysql_query(theConnection, theQuery);
        if (status) {
            std::cout << "Failed to perform the query.";
            return nullptr;
        }
        return mysql_store_result(theConnection);
    }
}