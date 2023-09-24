/// <summary>
/// Stores the constant parameters used in the library.
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#ifndef _GeometryParameters_HeaderFile
#define _GeometryParameters_HeaderFile

// Constants
namespace GeometryNamespace {
	enum DIMENSIONS
	{
		D2 = 2,
		D3 = 3
	};
	enum INTERSECTION1
	{
		Skew1 = -1,
		Intersects1 = 0,
		Coincides1 = 1
	};
	enum INTERSECTION2
	{
		Skew2 = -1,
		Intersects2 = 0,
		Includes2 = 1
	};
	const double TOLERANCE_GENERAL = 0.001;
	const double TOLERANCE_SENSITIVE = 0.000001;
}

#endif
