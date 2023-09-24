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
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
		Vector3D::Vector3D(const std::array<double, 3>& theLocalComponents)
		: VectorBase(DIMENSIONS::D3, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(const std::vector<double, std::allocator<double>>& theLocalComponents)
		: VectorBase(DIMENSIONS::D3, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// int* argument theNull is required to distinguish from the constructor: Vector3D(std::array<double, 3> theLocalComponents)
	/// Can have any value: null_ptr or 0 or 1 or 9999 or else
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Vector3D::Vector3D(
		const std::array<double, 3>& theAngles,
		bool&& anyValue)
		: VectorBase(DIMENSIONS::D3)
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// int* argument theNull is required to distinguish from the constructor: Vector3D(std::array<double, 3> theLocalComponents)
	/// Can have any value: null_ptr or 0 or 1 or 9999 or else
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Vector3D::Vector3D(
		const std::vector<double, std::allocator<double>>& theAngles,
		bool&& anyValue)
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
	Vector3D::Vector3D(ARGCOPY(PointBase) thePoint)
		: VectorBase(DIMENSIONS::D3, thePoint.getReferenceCoordSystem(), thePoint.getLocalCoords()) { }

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
		ARGCOPY(PointBase) thePoint1)
		: VectorBase(DIMENSIONS::D3, thePoint0, thePoint1) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalComponents)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector3D::Vector3D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalComponents)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// int* argument theNull is required to distinguish from the constructor:
	///	Vector3D(ARGCOPY(CoordSystem) theReferenceCoordSystem, std::array<double, 3> theLocalComponents)
	/// Can have any value: null_ptr or 0 or or 9999 or else
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector3D::Vector3D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theAngles,
		bool&& anyValue)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, std::array<double, 3>{{}})
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// int* argument theNull is required to distinguish from the constructor:
	///	Vector3D{CoordSystem8, theReferenceCoordSystem, const std::vector-double theLocalComponents}
	/// Can have any value: null_ptr or 0 or or 9999 or else
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector3D::Vector3D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theAngles,
		bool&& anyValue)
		: VectorBase(DIMENSIONS::D3, theReferenceCoordSystem, std::array<double, 3>{{}})
	{
		setLocalComponentsUsingAngles(theAngles);
	}

	/// <summary>
	/// Ctor
	/// Conversion from Vector2D to Vector3D
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector3D::Vector3D(ARGCOPY(Vector2D) theVector)
		: VectorBase(DIMENSIONS::D3)
	{
		c_dimensionCount = DIMENSIONS::D3;
		c_referenceCoordSystem = theVector.getReferenceCoordSystem();
		std::copy(
			std::begin(theVector.getLocalComponents()),
			std::end(theVector.getLocalComponents()),
			std::begin(c_localComponents));
	}

	/// <summary>
	/// Vector summation - Member function
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	template<typename T>
	auto Vector3D::operator+(T&& theVector) const->std::shared_ptr<Vector3D>
	{
		return add(std::forward<T>(theVector));
	}

	/// <summary>
	/// Vector subtruction - Member function
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	template<typename T>
	auto Vector3D::operator-(T&& theVector) const->std::shared_ptr<Vector3D>
	{
		return subtruct(std::forward<T>(theVector));
	}

	/// <summary>
	/// Getter - The unit vector in the same direction
	/// </summary>
	auto Vector3D::getUnitVector() const -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(
			c_referenceCoordSystem,
			getUnitVectorComponents());
	}

	/// <summary>
	/// Setter - Local components - By angles
	/// </summary>
	void Vector3D::setLocalComponentsUsingAngles(const std::array<double, 3>& theAngles)
	{
		auto localComponents = calculateComponentsFromAngles(theAngles);
		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// Setter - Local components - By angles
	/// </summary>
	void Vector3D::setLocalComponentsUsingAngles(const std::vector<double, std::allocator<double>>& theAngles)
	{
		if (theAngles.size() != DIMENSIONS::D3)
		{
			throw ArraySizeException();
		}

		std::array<double, 3> angles = { {} };
		std::copy(theAngles.begin(), theAngles.end(), angles.begin());
		auto localComponents = calculateComponentsFromAngles(angles);
		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// Create unit vector - Global- X
	/// </summary>
	auto Vector3D::createUnitVectorX() -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(std::array<double, 3>{{ 1., 0., 0. }});
	}

	/// <summary>
	/// Create unit vector - Global- Y
	/// </summary>
	auto Vector3D::createUnitVectorY() -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(std::array<double, 3>{{ 0., 1., 0. }});
	}

	/// <summary>
	/// Returns a vector resulted by a scalar multiplication
	/// Argument is a universal reference (rvalue reference is forwarded for performance)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	auto Vector3D::multiply(const double& theFactor) const -> std::shared_ptr<Vector3D>
	{
		if (GeometryMath::equals(theFactor, 0., getToleranceGeneral())) throw ZeroVectorException();

		auto finalCoords = GeometryMath::factorizeArray(c_localComponents, theFactor);
		return std::make_shared<Vector3D>(
			c_referenceCoordSystem,
			finalCoords);
	}
	/// <summary>
	/// Vector summation - Global function
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	template<typename T>
	auto operator+(ARGCOPY(Vector3D) theVector1, T&& theVector2)
		->std::shared_ptr<Vector3D>
	{
		return theVector1.add(std::forward<T>(theVector2));
	}

	/// <summary>
	/// Vector subtruction - Global function
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	template<typename T>
	auto operator-(ARGCOPY(Vector3D) theVector1, T&& theVector2)
		->std::shared_ptr<Vector3D>
	{
		return theVector1.subtruct(std::forward<T>(theVector2));
	}
}
