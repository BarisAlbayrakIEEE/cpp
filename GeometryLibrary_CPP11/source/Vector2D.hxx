/// <summary>
/// Vector2D defines a vector having z-component equal to zero wrt the reference CS.
/// 
/// INVARIANT:
///   1. A vector shall have a non-zero magnitude. Hence, at least one component shall be non-zero.
///   2. The local component in z-direction wrt the reference CS shall be zero.
/// 
/// Created for the tracebility and for the simplicity when switching to a local CS.
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

#ifndef _Vector2D_HeaderFile
#define _Vector2D_HeaderFile

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
	class Vector2D : public VectorBase
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

		// Private default ctor used for cloning the object
		Vector2D() = default;

	public:
		// ctor / dtor / operators
		explicit Vector2D(const std::array<double, 2>& theLocalComponents);
		explicit Vector2D(const std::array<double, 3>& theLocalComponents);
		explicit Vector2D(const std::vector<double, std::allocator<double>>& theLocalComponents);
		explicit Vector2D(const double& theAngle);
		explicit Vector2D(ARGCOPY(Point2D) thePoint);
		Vector2D(
			ARGCOPY(Point2D) thePoint0,
			ARGCOPY(Point2D) thePoint1);
		Vector2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 2>& theLocalComponents);
		Vector2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalComponents);
		Vector2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalComponents);
		Vector2D(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const double& theAngle);

		Vector2D(const Vector2D& rhs) = default;
		Vector2D& operator=(const Vector2D& rhs) = default;
		Vector2D(Vector2D&& rhs) noexcept = default;
		Vector2D& operator=(Vector2D&& rhs) noexcept = default;
		~Vector2D() final = default;

		template<typename T>
		auto operator+(T&& theVector) const->std::shared_ptr<T>;
		template<typename T>
		auto operator-(T&& theVector) const->std::shared_ptr<T>;

		// Methods
		auto getUnitVector() const -> std::shared_ptr<Vector2D>;
		double getSlope();
		double getAngle();
		auto createNormalVector();
		void setLocalComponentZ(const double& theLocalComponentZ) = delete;

		static auto createUnitVectorX();
		static auto createUnitVectorY();
		template<typename T>
		auto add(T&& theVector) const
			-> std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>>;
		template<typename T>
		auto subtruct(T&& theVector) const
			-> std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector)>>>;
		auto multiply(const double& theFactor) const->std::shared_ptr<Vector2D>;
	};

	template<typename T>
	auto operator+(ARGCOPY(Vector2D) theVector1, T&& theVector2)
		->std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector2)>>>;
	template<typename T>
	auto operator-(ARGCOPY(Vector2D) theVector1, T&& theVector2)
		->std::shared_ptr<typename std::remove_const_t<std::remove_reference_t<decltype(theVector2)>>>;
}

#endif
