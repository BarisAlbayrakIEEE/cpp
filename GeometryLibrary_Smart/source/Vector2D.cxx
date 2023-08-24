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

	/// <summary>
	/// Ctor
	/// Vector by components
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const arrayS2& theLocalComponents)
		: VectorBase(DIMENSIONS::D2)
	{
		arrayS3 localComponents{ GeometryMath::convertArrayS2ToS3(theLocalComponents) };
		inspectLocalComponents(localComponents);
		c_localComponents = localComponents;
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const arrayS3& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, *GlobalCoordSystem::getGlobalCoordSystem(), theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const vectorInput1D& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, *GlobalCoordSystem::getGlobalCoordSystem(), theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Vector2D::Vector2D(const double& theAngle)
		: VectorBase(DIMENSIONS::D2, *GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{ 1., theAngle, 0. }) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the reference CS of the point
	/// Components are defined by a point
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(ARGCOPY(Point2D) thePoint)
		: VectorBase(DIMENSIONS::D2, thePoint.getReferenceCoordSystem(), thePoint.getLocalCoords()) { }

	/// <summary>
	/// Ctor
	/// Vector by two points
	/// Reference CS is the reference CS of the points
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		ARGCOPY(Point2D) thePoint0,
		ARGCOPY(Point2D) thePoint1)
		: VectorBase(DIMENSIONS::D2, thePoint0, thePoint1) { }

	/// <summary>
	/// Ctor
	/// Vector by components
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		CoordSystem& theReferenceCoordSystem,
		const arrayS2& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem)
	{
		arrayS3 localComponents{ GeometryMath::convertArrayS2ToS3(theLocalComponents) };
		inspectLocalComponents(localComponents);
		c_localComponents = localComponents;
	}

	/// <summary>
	/// Ctor
	/// Vector by components
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		CoordSystem& theReferenceCoordSystem,
		const arrayS3& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Vector by components
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		CoordSystem& theReferenceCoordSystem,
		const vectorInput1D& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Components are defined by an angle
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector2D::Vector2D(
		CoordSystem& theReferenceCoordSystem,
		const double& theAngle)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem, arrayS3{ 1., theAngle, 0. }) { }

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector2D::Vector2D(const Vector2D& rhs)
		: VectorBase(DIMENSIONS::D2)
	{
		VectorBase::copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector2D& Vector2D::operator=(const Vector2D& rhs)
	{
		if (&rhs == this) return *this;

		VectorBase::copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector2D::Vector2D(Vector2D&& rhs) noexcept
		: VectorBase(DIMENSIONS::D2)
	{
		VectorBase::copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector2D& Vector2D::operator=(Vector2D&& rhs) noexcept
	{
		if (&rhs == this) return *this;

		VectorBase::copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	Vector2D::~Vector2D() {
		Destroy();
	}

	/// <summary>
	/// Getter - The unit vector in the same direction
	/// </summary>
	Vector2D Vector2D::getUnitVector() const {
		return Vector2D(*c_referenceCoordSystem, getUnitVectorComponents());
	}

	double Vector2D::getSlope() {
		return calculateSlope(c_localComponents[0], c_localComponents[1]);
	}

	double Vector2D::getAngle() {
		return calculateAngle(c_localComponents[0], c_localComponents[1]);
	}

	Vector2D Vector2D::createNormalVector() {
		return Vector2D(*c_referenceCoordSystem, arrayS3{ 1., -1. / getAngle(), 0. });
	}
}
