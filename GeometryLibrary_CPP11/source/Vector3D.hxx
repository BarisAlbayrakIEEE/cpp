/// <summary>
/// Vector3D defines a vector wrt the reference CS.
/// 
/// INVARIANT:
///   A vector shall have a non-zero magnitude.
///   Hence, at least one component shall be non-zero.
/// 
/// Created for the tracebility and for the simplicity when switching between CSs.
/// 
/// Uses universal references in add and subtruct functions
/// which allow rvalue references can be utilized together with the lvalue references
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Vector3D_HeaderFile
#define _Vector3D_HeaderFile

#ifndef _VectorBase_HeaderFile
#include "VectorBase.hxx"
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
	class Vector3D : public VectorBase
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

		// Private default ctor used for cloning the object
		Vector3D() = default;

	public:
		// ctor / dtor / operators
		explicit Vector3D(const std::array<double, 3>& theLocalComponents);
		explicit Vector3D(const std::vector<double, std::allocator<double>>& theLocalComponents);
		Vector3D(const std::array<double, 3>& theAngles, bool&& anyValue);
		Vector3D(const std::vector<double, std::allocator<double>>& theAngles, bool&& anyValue);
		explicit Vector3D(ARGCOPY(PointBase) thePoint);
		Vector3D(
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1);
		Vector3D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalComponents);
		Vector3D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalComponents);
		Vector3D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theAngles,
			bool&& anyValue);
		Vector3D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theAngles,
			bool&& anyValue);
		explicit Vector3D(ARGCOPY(Vector2D) theVector);

		Vector3D(const Vector3D& rhs) = default;
		Vector3D& operator=(const Vector3D& rhs) = default;
		Vector3D(Vector3D&& rhs) noexcept = default;
		Vector3D& operator=(Vector3D&& rhs) noexcept = default;
		~Vector3D() final = default;

		template<typename T>
		auto operator+(T&& theVector) const->std::shared_ptr<Vector3D>;
		template<typename T>
		auto operator-(T&& theVector) const->std::shared_ptr<Vector3D>;

		// Methods
		auto getUnitVector() const -> std::shared_ptr<Vector3D>;
		void setLocalComponentsUsingAngles(const std::array<double, 3>& theAngles);
		void setLocalComponentsUsingAngles(const std::vector<double, std::allocator<double>>& theAngles);

		static auto createUnitVectorX() -> std::shared_ptr<Vector3D>;
		static auto createUnitVectorY() -> std::shared_ptr<Vector3D>;

		template<typename T>
		inline auto add(T&& theVector) const
			->std::shared_ptr<Vector3D>
		{
			// Inspect if a Vector object
			static_assert(
				std::is_base_of<VectorBase, typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>>::value,
				"T must inherit from VectorBase");

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
			return std::make_shared<Vector3D>(
				c_referenceCoordSystem,
				localCoords);
		};

		template<typename T>
		inline auto subtruct(T&& theVector) const
			->std::shared_ptr<Vector3D>
		{
			// Inspect if a Vector object
			static_assert(
				std::is_base_of<VectorBase, typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>>::value,
				"T must inherit from VectorBase");

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
			return std::make_shared<Vector3D>(
				c_referenceCoordSystem,
				localCoords);
		};

		auto multiply(const double& theFactor) const -> std::shared_ptr<Vector3D>;
	};

	template<typename T>
	auto operator+(ARGCOPY(Vector3D) theVector1, T&& theVector2)
		->std::shared_ptr<Vector3D>;
	template<typename T>
	auto operator-(ARGCOPY(Vector3D) theVector1, T&& theVector2)
		->std::shared_ptr<Vector3D>;
}

#endif
