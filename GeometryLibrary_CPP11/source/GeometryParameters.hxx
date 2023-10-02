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
	class GeometryParameters {
		static double TOLERANCE_GENERAL;
		static double TOLERANCE_SENSITIVE;

	public:
		static double getToleranceGeneral()
		{
			return TOLERANCE_GENERAL;
		};
		static double getToleranceSensitive()
		{
			return TOLERANCE_SENSITIVE;
		};
		static void setToleranceGeneral(const double& theToleranceGeneral)
		{
			TOLERANCE_GENERAL = theToleranceGeneral;
		};
		static void setToleranceSensitive(const double& theToleranceSensitive)
		{
			TOLERANCE_SENSITIVE = theToleranceSensitive;
		};

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
	};
}

#endif
