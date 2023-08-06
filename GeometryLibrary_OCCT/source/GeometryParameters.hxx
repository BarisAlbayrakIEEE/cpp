/// <summary>
/// Stores the constant parameters used in the library.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#ifndef _GeometryParameters_HeaderFile
#define _GeometryParameters_HeaderFile

#define _USE_MATH_DEFINES
#include <math.h>

// Constants
namespace GeometryNamespace {
	const enum OBJECT_TYPES
	{
		CoordSystem_,
		Point2D_,
		Point3D_,
		Vector2D_,
		Vector3D_,
		Axis_,
		Line_,
		Circle_,
		Plane_
	};
	const enum DIMENSIONS
	{
		D2 = 2,
		D3 = 3
	};
	const enum INTERSECTION1
	{
		Skew1 = -1,
		Intersects1 = 0,
		Coincides1 = 1
	};
	const enum INTERSECTION2
	{
		Skew2 = -1,
		Intersects2 = 0,
		Includes2 = 1
	};
	const double TOLERANCE_GENERAL = 0.001;
	const double TOLERANCE_SENSITIVE = 0.000001;
	const double INFINITE_VALUE = 1.E20;
}

#endif
