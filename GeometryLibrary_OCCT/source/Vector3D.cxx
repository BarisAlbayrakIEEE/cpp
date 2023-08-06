/// <summary>
/// Ctor
/// Reference CS is the global CS
/// baris.albayrak.ieee@gmail.com
/// </summary>

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
	IMPLEMENT_STANDARD_RTTIEXT(Vector3D, VectorBase)

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
		Vector3D::Vector3D(const arrayS3& theLocalComponents) throw (ZeroVectorException)
		: VectorBase(DIMENSIONS::D3, GlobalCoordSystem::getGlobalCoordSystem(), theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(const vectorInput1D& theLocalComponents) throw (ZeroVectorException)
		: VectorBase(DIMENSIONS::D3, GlobalCoordSystem::getGlobalCoordSystem(), theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// int argument theNull is required to distinguish from the constructor: Vector3D(arrayS3 theLocalComponents)
	/// Can have any value: Nane or 0 or 1or 9999 or else
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Vector3D::Vector3D(const arrayS3& theAngles, int* theNull)
		: VectorBase(DIMENSIONS::D3)
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// int argument theNull is required to distinguish from the constructor: Vector3D(arrayS3 theLocalComponents)
	/// Can have any value: Nane or 0 or 1or 9999 or else
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Vector3D::Vector3D(const vectorInput1D& theAngles, int* theNull)
		: VectorBase(DIMENSIONS::D3)
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the reference CS of the point
	/// Components are defined by a point
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(ARGCOPY(PointBase) thePoint) throw (NullptrException, ZeroVectorException)
		: VectorBase(DIMENSIONS::D3, thePoint->getReferenceCoordSystem(), thePoint->getLocalCoords()) { }

	/// <summary>
	/// Ctor
	/// Vector by two points
	/// Reference CS is the reference CS of the points
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(
		ARGCOPY(PointBase) thePoint0,
		ARGCOPY(PointBase) thePoint1) throw (NullptrException, CoordSystemMismatchException, ZeroVectorException)
		: VectorBase(DIMENSIONS::D3, thePoint0, thePoint1) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theLocalComponents) throw (NullptrException, ZeroVectorException)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theLocalComponents) throw (NullptrException, ZeroVectorException)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// int argument theNull is required to distinguish from the constructor:
	///	Vector3D(CoordSystem& theReferenceCoordSystem, arrayS3 theLocalComponents)
	/// Can have any value: Nane or 0 or 1or 9999 or else
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector3D::Vector3D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theAngles, int* theNull) throw (NullptrException)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, arrayS3{})
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// int argument theNull is required to distinguish from the constructor:
	///	Vector3D{CoordSystem8, theReferenceCoordSystem, const std::vector-double theLocalComponents}
	/// Can have any value: None or 0 or 1or 9999 or else
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector3D::Vector3D(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theAngles, int* theNull) throw (NullptrException)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, arrayS3{})
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// Ctor
	/// Conversion from Vector2D to Vector3D
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector3D::Vector3D(ARGCOPY(Vector2D) theVector) throw (NullptrException)
		: VectorBase(DIMENSIONS::D3)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		c_dimensionCount = DIMENSIONS::D3;
		c_referenceCoordSystem = theVector->getReferenceCoordSystem();
		std::copy(std::begin(theVector->getLocalComponents()), std::end(theVector->getLocalComponents()), std::begin(c_localComponents));
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector3D::Vector3D(const Vector3D& rhs)
		: VectorBase(DIMENSIONS::D3)
	{
		VectorBase::copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector3D& Vector3D::operator=(const Vector3D& rhs)
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
	Vector3D::Vector3D(Vector3D&& rhs) noexcept
		: VectorBase(DIMENSIONS::D3)
	{
		VectorBase::copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Vector3D& Vector3D::operator=(Vector3D&& rhs) noexcept
	{
		if (&rhs == this) return *this;

		VectorBase::copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// Use OCCT approach
	/// </summary>
	Vector3D::~Vector3D() {
		Destroy();
	}

	/// <summary>
	/// Getter - The unit vector in the same direction
	/// </summary>
	OUTVAL(Vector3D) Vector3D::getUnitVector() const {
		return new Vector3D(c_referenceCoordSystem, getUnitVectorComponents());
	}

	/// <summary>
	/// Setter - Local components - By angles
	/// </summary>
	void Vector3D::setLocalComponentsUsingAngles(const arrayS3& theAngles)
	{
		arrayS3 localComponents = calculateComponentsFromAngles(theAngles);
		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// Setter - Local components - By angles
	/// </summary>
	void Vector3D::setLocalComponentsUsingAngles(const vectorInput1D& theAngles) throw (ArraySizeException)
	{
		if (theAngles.size() != DIMENSIONS::D3)
		{
			throw ArraySizeException();
		}

		arrayS3 angles{ GeometryMath::convertVectorToArray1DS3(theAngles) };
		arrayS3 localComponents = calculateComponentsFromAngles(angles);
		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}
}
