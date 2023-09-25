// baris.albayrak.ieee@gmail.com

#include "Macros.h"
#include "GeometryObject.hxx"
#include "GeometryParameters.hxx"
#include "GeometryMath.hxx"
#include "GeometryException.hxx"
#include "ReferenceObject.hxx"
#include "CoordSystem.hxx"
#include "GlobalCoordSystem.hxx"
#include "PointBase.hxx"
#include "Point2D.hxx"
#include "Point3D.hxx"
#include "VectorBase.hxx"
#include "Vector2D.hxx"
#include "Vector3D.hxx"
#include "Axis.hxx"
#include "Line.hxx"
#include "Circle.hxx"
#include "Plane.hxx"

namespace GeometryNamespace {

	/// <summary>
	/// The default constructor
	/// Reference CS is the global CS
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D()
		: PointBase(DIMENSIONS::D2, std::array<double, 3>{{}}) { }

	/// <summary>
	/// Ctor
	/// Point by coords
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D(const std::array<double, 2>& theLocalComponents)
		: PointBase(DIMENSIONS::D2)
	{
		std::copy(theLocalComponents.begin(), theLocalComponents.end(), c_localCoords.begin());
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D(const std::array<double, 3>& theLocalCoords)
		: PointBase(DIMENSIONS::D2, theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D(const std::vector<double, std::allocator<double>>& theLocalCoords)
		: PointBase(DIMENSIONS::D2, theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem, std::array<double, 3>{{}}) { }

	/// <summary>
	/// Ctor
	/// Point by coords
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 2>& theLocalComponents)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem)
	{
		std::copy(theLocalComponents.begin(), theLocalComponents.end(), c_localCoords.begin());
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalCoords)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalCoords)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalCoords) { }
}
