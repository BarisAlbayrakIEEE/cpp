/// <summary>
/// Point3D defines a point wrt the reference CS.
/// 
/// NO INVARIANT
/// 
/// Created for the tracebility and for the simplicity when switching between CSs.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
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
		Point3D();
		explicit Point3D(const std::array<double, 3>& theLocalCoords);
		explicit Point3D(const std::vector<double, std::allocator<double>>& theLocalCoords);
		explicit Point3D(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
		Point3D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalCoords);
		Point3D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalCoords);
		explicit Point3D(ARGCOPY(Point2D) thePoint);

		Point3D(const Point3D& rhs) = default;
		Point3D& operator=(const Point3D& rhs) = default;
		Point3D(Point3D&& rhs) noexcept = default;
		Point3D& operator=(Point3D&& rhs) noexcept = default;
		~Point3D() final = default;
	};
}

#endif
