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
	IMPLEMENT_STANDARD_RTTIEXT(Point3D, PointBase)

	/// <summary>
	/// The default constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point3D::Point3D()
		: PointBase(DIMENSIONS::D3, GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{}) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point3D::Point3D(const arrayS3& theLocalCoords)
		: PointBase(DIMENSIONS::D3, GlobalCoordSystem::getGlobalCoordSystem(), theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Point3D::Point3D(const vectorInput1D& theLocalCoords)
		: PointBase(DIMENSIONS::D3, GlobalCoordSystem::getGlobalCoordSystem(), theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(ARGCOPY(CoordSystem) theReferenceCoordSystem)
		: PointBase(DIMENSIONS::D3, theReferenceCoordSystem, arrayS3{}) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theLocalCoords)
		: PointBase(DIMENSIONS::D3, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theLocalCoords)
		: PointBase(DIMENSIONS::D3, theReferenceCoordSystem, theLocalCoords) { }

	/// <summary>
	/// Ctor
	/// Conversion from Point2D to Point3D
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Point3D::Point3D(ARGCOPY(Point2D) thePoint)
		: PointBase(DIMENSIONS::D3)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		c_dimensionCount = DIMENSIONS::D3;
		c_referenceCoordSystem = thePoint->getReferenceCoordSystem();
		std::copy(std::begin(thePoint->getLocalCoords()), std::end(thePoint->getLocalCoords()), std::begin(c_localCoords));
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point3D::Point3D(const Point3D& rhs)
		: PointBase(rhs)
	{
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point3D& Point3D::operator=(const Point3D& rhs)
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
	Point3D::Point3D(Point3D&& rhs) noexcept
		: PointBase(rhs)
	{
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Point3D& Point3D::operator=(Point3D&& rhs) noexcept
	{
		if (&rhs == this) return *this;

		PointBase::copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// Use OCCT approach
	/// </summary>
	Point3D::~Point3D() {
		Destroy();
	}
}
