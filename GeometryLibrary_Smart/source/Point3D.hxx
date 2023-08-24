/// <summary>
/// Point3D defines a point wrt the reference CS.
/// 
/// NO INVARIANT
/// 
/// Created for the tracebility and for the simplicity when switching between CSs.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Point3D_HeaderFile
#define _Point3D_HeaderFile

#ifndef _PointBase_HeaderFile
#include "PointBase.hxx"
#endif
#ifndef _GeometryParameters_HeaderFile
#include "GeometryParameters.hxx"
#endif
#ifndef _GeometryMath_HeaderFile
#include "GeometryMath.hxx"
#endif
#ifndef _GeometryException_HeaderFile
#include "GeometryException.hxx"
#endif
#ifndef _Macros_HeaderFile
#include "Macros.h"
#endif

namespace GeometryNamespace {
	class GeometryObject;
	class ReferenceObject;
	class CoordSystem;
	class GlobalCoordSystem;
	class PointBase;
	class Point2D;
	class Point3D;
	class VectorBase;
	class Vector2D;
	class Vector3D;
	class Axis;
	class Line;
	class Circle;
	class Plane;

	class Point3D : public PointBase
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

	public:
		// ctor / dtor / operators
		Point3D(const Point3D& rhs);
		Point3D& operator=(const Point3D& rhs);
		Point3D(Point3D&& rhs) noexcept;
		Point3D& operator=(Point3D&& rhs) noexcept;
		~Point3D();

		Point3D();
		Point3D(const arrayS3& theLocalCoords);
		Point3D(const vectorInput1D& theLocalCoords);
		Point3D(CoordSystem& theReferenceCoordSystem);
		Point3D(
			CoordSystem& theReferenceCoordSystem,
			const arrayS3& theLocalCoords);
		Point3D(
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords);
		Point3D(ARGCOPY(Point2D) thePoint);
	};
}

#endif
