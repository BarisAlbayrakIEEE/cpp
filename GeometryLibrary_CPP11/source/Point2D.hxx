/// <summary>
/// Point2D defines a point having z-coordinate equal to zero wrt the reference CS.
/// 
/// INVARIANT:
///   The local coordinate in z-direction wrt the reference CS shall be zero.
/// 
/// Created for the tracebility and for the simplicity when switching to a local CS.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Point2D_HeaderFile
#define _Point2D_HeaderFile

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
	class Point2D : public PointBase
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

	public:
		// ctor / dtor / operators
		Point2D();
		explicit Point2D(const std::array<double, 2>& theLocalCoords);
		explicit Point2D(const std::array<double, 3>& theLocalCoords);
		explicit Point2D(const std::vector<double, std::allocator<double>>& theLocalCoords);
		explicit Point2D(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
		Point2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 2>& theLocalCoords);
		Point2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalCoords);
		Point2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalCoords);

		Point2D(const Point2D& rhs) = default;
		Point2D& operator=(const Point2D& rhs) = default;
		Point2D(Point2D&& rhs) noexcept = default;
		Point2D& operator=(Point2D&& rhs) noexcept = default;
		~Point2D() final = default;

		// Methods
		void setLocalCoordZ(double theLocalCoordZ) = delete;
	};
}

#endif
