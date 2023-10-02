// baris.albayrak.ieee@gmail.com

#include "Macros.h"
#include "GeometryObject.hxx"
#include "GeometryMath.hxx"
#include "GeometryParameters.hxx"
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
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point3D::Point3D()
		: PointBase(GeometryParameters::DIMENSIONS::D3, std::array<double, 3>{{}}) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point3D::Point3D(const std::array<double, 3>& theLocalCoords)
		: PointBase(GeometryParameters::DIMENSIONS::D3, theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point3D::Point3D(const std::vector<double, std::allocator<double>>& theLocalCoords)
		: PointBase(GeometryParameters::DIMENSIONS::D3, theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		: PointBase(GeometryParameters::DIMENSIONS::D3, theReferenceCoordSystem, std::array<double, 3>{{}}) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalCoords)
		: PointBase(GeometryParameters::DIMENSIONS::D3, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalCoords)
		: PointBase(GeometryParameters::DIMENSIONS::D3, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Conversion from Point2D to Point3D
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(ARGCOPY(Point2D) thePoint)
		: PointBase(GeometryParameters::DIMENSIONS::D3)
	{
		c_dimensionCount = GeometryParameters::DIMENSIONS::D3;
		c_referenceCoordSystem = thePoint.getReferenceCoordSystem();
		std::copy(std::begin(thePoint.getLocalCoords()), std::end(thePoint.getLocalCoords()), std::begin(c_localCoords));
	}
}
