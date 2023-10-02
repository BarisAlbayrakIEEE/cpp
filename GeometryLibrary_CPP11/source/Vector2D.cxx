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
	/// Ctor
	/// Vector by components
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const std::array<double, 2>& theLocalComponents)
		: VectorBase(GeometryParameters::DIMENSIONS::D2)
	{
		std::array<double, 3> localComponents = { {} };
		std::copy(theLocalComponents.cbegin(), theLocalComponents.cend(), localComponents.begin());
		inspectLocalComponents(localComponents.cbegin(), localComponents.cend());
		setLocalComponents(localComponents);
	}

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const std::array<double, 3>& theLocalComponents)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(const std::vector<double, std::allocator<double>>& theLocalComponents)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// Components are defined by an angle
	/// </summary>
	Vector2D::Vector2D(const double& theAngle)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, std::array<double, 3>{ 1., theAngle, 0. }) { }

	/// <summary>
	/// Ctor
	/// Reference CS is the reference CS of the point
	/// Components are defined by a point
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(ARGCOPY(Point2D) thePoint)
		: VectorBase(
			GeometryParameters::DIMENSIONS::D2,
			thePoint.getReferenceCoordSystem(),
			thePoint.getLocalCoords()) { }

	/// <summary>
	/// Ctor
	/// Vector by two points
	/// Reference CS is the reference CS of the points
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		ARGCOPY(Point2D) thePoint0,
		ARGCOPY(Point2D) thePoint1)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, thePoint0, thePoint1) { }

	/// <summary>
	/// Ctor
	/// Vector by components
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 2>& theLocalComponents)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, theReferenceCoordSystem)
	{
		std::array<double, 3> localComponents = { {} };
		std::copy(theLocalComponents.begin(), theLocalComponents.end(), localComponents.begin());
		inspectLocalComponents(localComponents.cbegin(), localComponents.cend());
		setLocalComponents(localComponents);
	}

	/// <summary>
	/// Ctor
	/// Vector by components
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalComponents)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Vector by components
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	Vector2D::Vector2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalComponents)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, theReferenceCoordSystem, theLocalComponents) { }

	/// <summary>
	/// Ctor
	/// Components are defined by an angle
	/// </summary>
	/// <exception> NullptrException </exception>
	Vector2D::Vector2D(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const double& theAngle)
		: VectorBase(GeometryParameters::DIMENSIONS::D2, theReferenceCoordSystem, std::array<double, 3>{ 1., theAngle, 0. }) { }

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
		std::array<double, 3> localCoords{ {} };
		std::transform(
			c_localComponents.cbegin(),
			c_localComponents.cend(),
			getVectorWithMyCoordSystem(std::forward<T>(theVector))->getLocalComponents().cbegin(),
			localCoords.begin(),
			std::plus<double>());
		if (
			std::all_of(
				localCoords.cbegin(),
				localCoords.cend(),
				[](double i) { return GeometryMath::zero_g(i); })) {
			throw ZeroVectorException();
		}
		if (theVector.is2D() && !localCoords[2].zero()) {
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
		std::array<double, 3> localCoords{ {} };
		std::transform(
			c_localComponents.cbegin(),
			c_localComponents.cend(),
			getVectorWithMyCoordSystem(std::forward<T>(theVector))->getLocalComponents().cbegin(),
			localCoords.begin(),
			std::minus<double>());
		if (
			std::all_of(
				localCoords.cbegin(),
				localCoords.cend(),
				[](double i) { return GeometryMath::zero_g(i); })) {
			throw ZeroVectorException();
		}
		if (theVector.is2D() && !localCoords[2].zero()) {
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
		if (GeometryMath::zero_g(theFactor)) throw ZeroVectorException();

		std::array<double, 3> finalCoords = { {} };
		std::transform(
			c_localComponents.begin(),
			c_localComponents.end(),
			finalCoords.begin(),
			[theFactor](double i) { return i * theFactor; });

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
