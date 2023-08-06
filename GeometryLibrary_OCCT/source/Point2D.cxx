// baris.albayrak.ieee@gmail.com

#include "GeometryObject.hxx"
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
	IMPLEMENT_STANDARD_RTTIEXT(Point2D, PointBase)

	/// <summary>
	/// The default constructor
	/// Reference CS is the global CS
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D()
		: PointBase(DIMENSIONS::D2, GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{}) { }

	/// <summary>
	/// Ctor
	/// Point by coords
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D(const arrayS2& theLocalComponents)
		: PointBase(DIMENSIONS::D2)
	{
		c_localCoords = GeometryMath::convertArrayS2ToS3(theLocalComponents);
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D(const arrayS3& theLocalCoords)
		: PointBase(DIMENSIONS::D2, GlobalCoordSystem::getGlobalCoordSystem(), theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point2D::Point2D(const vectorInput1D& theLocalCoords)
		: PointBase(DIMENSIONS::D2, GlobalCoordSystem::getGlobalCoordSystem(), theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(ARGCOPY(CoordSystem) theReferenceCoordSystem) throw (NullptrException)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem, arrayS3{}) { }

	/// <summary>
	/// Ctor
	/// Point by coords
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS2& theLocalComponents) throw (NullptrException)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem)
	{
		c_localCoords = GeometryMath::convertArrayS2ToS3(theLocalComponents);
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theLocalCoords) throw (NullptrException)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point2D::Point2D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theLocalCoords) throw (NullptrException)
		: PointBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point2D::Point2D(const Point2D& rhs)
		: PointBase(rhs)
	{
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point2D& Point2D::operator=(const Point2D& rhs)
	{
		if (&rhs == this) return *this;

		PointBase::copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point2D::Point2D(Point2D&& rhs) noexcept
		: PointBase(rhs)
	{
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point2D& Point2D::operator=(Point2D&& rhs) noexcept
	{
		if (&rhs == this) return *this;

		PointBase::copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// Use OCCT approach
	/// </summary>
	Point2D::~Point2D() {
		Destroy();
	}
}
