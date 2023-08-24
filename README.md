# cpp
C++ projects

>
>
>
>
GeometryLibrary_OCCT:
Simulates a simple geometry library (point, circle etc.) for which the objetc lifetime is managed by Standard_Handles defined in OpenCascade (OCCT) library.
This library is quite an old library (before 2011) for which the smart pointers (C++11) are not published yet.
Fr details of the project, see the docstring in GeometryObject.hxx which is the header file for the base object for the project.
Current release contains RAII issues (memory leaks) due to the complexity of OCCT Standard_Handle.
The problem will be solved in a future release.

>
>
>
>
GeometryLibrary_Smart:
GeometryLibrary_OCCT is re-implemented to use the shared pointers to simulate a strong ownership between the objects.
For example, a line object is defined by 2 point objects and makes an association relation with the axis class.
Hence, the resource acquired for the two points and the associated axis should be released when the line (and other objects sharing the ownership) is out of scope.
With Standard_Handles of OCCT its hard to implement this relation when there is a hierarchy between the classes
as some of the compilation errors about OCCT handles (e.g. a call to a deleted copy assignment of a base class) can be missed by compilers (e.g. gcc).
Hence, RAII problems occurs in runtime.
Actually, the available revision of GeometryLibrary_OCCT contains such memory leaks currently.
A future release will cover this issue.

>
>
>
MYSQL_Library:
This is a simple library containinf a few methods such as creating a connection and performing queries.
