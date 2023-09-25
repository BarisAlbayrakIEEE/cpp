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
	/// Vector by components
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const std::array<double, 2>& theLocalComponents)
		: VectorBase(DIMENSIONS::D2)
	{
		std::array<double, 3> localComponents = { {} };
		std::copy(theLocalComponents.begin(), theLocalComponents.end(), localComponents.begin());
		inspectLocalComponents(localComponents);
		c_localComponents = localComponents;
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const std::array<double, 3>& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const std::vector<double, std::allocator<double>>& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	Vector2D::Vector2D(const double& theAngle)
		: VectorBase(DIMENSIONS::D2, std::array<double, 3>{ 1., theAngle, 0. }) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the reference CS of the point
	/// Components are defined by a point
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(ARGCOPY(Point2D) thePoint)
		: VectorBase(
			DIMENSIONS::D2,
			thePoint.getReferenceCoordSystem(),
			thePoint.getLocalCoords()) { }

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
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 2>& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem)
	{
		std::array<double, 3> localComponents = { {} };
		std::copy(theLocalComponents.begin(), theLocalComponents.end(), localComponents.begin());
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
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Vector by components
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalComponents)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Components are defined by an angle
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector2D::Vector2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const double& theAngle)
		: VectorBase(DIMENSIONS::D2, theReferenceCoordSystem, std::array<double, 3>{ 1., theAngle, 0. }) { }

	/// <summary>
	/// Vector summation - Member function
	/// Argument is a universal reference (rvalue reference is forwarded for performance)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	template<typename T>
	auto Vector2D::operator+(T&& theVector) const->std::shared_ptr<T>
	{
		return add(std::forward<T>(theVector));
	}

	/// <summary>
	/// Vector subtruction - Member function
	/// Argument is a universal reference (rvalue reference is forwarded for performance)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	template<typename T>
	auto Vector2D::operator-(T&& theVector) const->std::shared_ptr<T>
	{
		return subtruct(std::forward<T>(theVector));
	}


	/// <summary>
	/// Getter - The unit vector in the same direction
	/// </summary>
	auto Vector2D::getUnitVector() const -> std::shared_ptr<Vector2D>
	{
		return std::make_shared<Vector2D>(
			c_referenceCoordSystem,
			getUnitVectorComponents());
	}

	double Vector2D::getSlope() {
		return calculateSlope(c_localComponents[0], c_localComponents[1]);
	}

	double Vector2D::getAngle() {
		return calculateAngle(c_localComponents[0], c_localComponents[1]);
	}

	auto Vector2D::createNormalVector() {
		return std::make_shared<Vector2D>(
			c_referenceCoordSystem,
			std::array<double, 3>{ 1., -1. / getAngle(), 0. });
	}

	/// <summary>
	/// Create unit vector - Global- X
	/// </summary>
	auto Vector2D::createUnitVectorX() {
		return std::make_shared<Vector2D>(std::array<double, 3>{ 1., 0., 0. });
	}

	/// <summary>
	/// Create unit vector - Global- Y
	/// </summary>
	auto Vector2D::createUnitVectorY()
	{
		return std::make_shared<Vector2D>(std::array<double, 3>{ 0., 1., 0. });
	}

	/// <summary>
	/// Returns a vector resulted by the sum with the input vector (vector summation)
	/// Returns Vector2D if the input vector is Vector2D
	/// Otherwise, return Vector3D.
	/// Argument is a universal reference (rvalue reference is forwarded for performance)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	template<typename T>
	auto Vector2D::add(T&& theVector) const
		-> std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>>
	{
		using _typeName = typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>;

		// Inspect if a Vector object
		static_assert(std::is_base_of<VectorBase, _typeName>::value, "T must inherit from VectorBase");

		// Get the item in my reference CS
		std::array<double, 3> localCoords{
			GeometryMath::sumArrays(
				c_localComponents,
				getVectorWithMyCoordSystem(std::forward<T>(theVector))->getLocalComponents()
			)
		};
		if (GeometryMath::equalsZero(localCoords, getToleranceGeneral())) {
			throw ZeroVectorException();
		}
		if (theVector.is2D() && !GeometryMath::equals(localCoords[2], 0., getToleranceGeneral())) {
			throw CoordSystemMismatchException();
		}
		return std::make_shared<_typeName>(
			c_referenceCoordSystem,
			localCoords);
	}

	/// <summary>
	/// Returns a vector resulted by the sum with the reverse of the input vector (vector subtruction)
	/// Returns Vector2D if the input vector is Vector2D
	/// Otherwise, return Vector3D.
	/// Argument is a universal reference (rvalue reference is forwarded for performance)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	template<typename T>
	auto Vector2D::subtruct(T&& theVector) const
		-> std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>>
	{
		using _typeName = typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>;

		// Inspect if a Vector object
		static_assert(std::is_base_of<VectorBase, _typeName>::value, "T must inherit from VectorBase");

		// Get the item in my reference CS
		std::array<double, 3> localCoords{
			GeometryMath::subtructArrays(
				c_localComponents,
				getVectorWithMyCoordSystem(std::forward<T>(theVector))->getLocalComponents()
			)
		};
		if (GeometryMath::equalsZero(localCoords, getToleranceGeneral())) {
			throw ZeroVectorException();
		}
		if (theVector.is2D() && !GeometryMath::equals(localCoords[2], 0., getToleranceGeneral())) {
			throw CoordSystemMismatchException();
		}
		return std::make_shared<_typeName>(
			c_referenceCoordSystem,
			localCoords);
	}

	/// <summary>
	/// Returns a vector resulted by a scalar multiplication
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	auto Vector2D::multiply(const double& theFactor) const -> std::shared_ptr<Vector2D>
	{
		if (GeometryMath::equals(theFactor, 0., getToleranceGeneral())) throw ZeroVectorException();

		auto finalCoords = GeometryMath::factorizeArray(c_localComponents, theFactor);
		return std::make_shared<Vector2D>(
			c_referenceCoordSystem,
			finalCoords);
	}
	/// <summary>
	/// Vector summation - Global function
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	template<typename T>
	auto operator+(ARGCOPY(Vector2D) theVector1, T&& theVector2)
		->std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector2)>>>
	{
		return theVector1.add(std::forward<T>(theVector2));
	}

	/// <summary>
	/// Vector subtruction - Global function
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	template<typename T>
	auto operator-(ARGCOPY(Vector2D) theVector1, T&& theVector2)
		->std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector2)>>>
	{
		return theVector1.subtruct(std::forward<T>(theVector2));
	}
}
