/// <summary>
/// VectorBase defines a vector and its members.
/// 
/// INVARIANT:
///   A vector shall have a non-zero magnitude.
///   Hence, at least one component shall be non-zero.
/// 
/// Vector2D and Vector3D are the child classes.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4250)
#pragma warning(disable : 4290)

#ifndef _VectorBase_HeaderFile
#define _VectorBase_HeaderFile

#ifndef _ReferenceObject_HeaderFile
#include "ReferenceObject.hxx"
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

namespace GeometryNamespace{
	class VectorBase : public ReferenceObject
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

		using ReferenceObject::is2D; // To disable 4250 warning
		using ReferenceObject::is3D; // To disable 4250 warning

		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		std::array<double, 3> c_localComponents = { {} };

	protected:
		// Private default ctor used for cloning the object
		VectorBase() = default;

		// ctor / dtor / operators
		explicit VectorBase(const int theDimensionCount);
		VectorBase(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
		VectorBase(
			const int theDimensionCount,
			const std::array<double, 3>& theLocalComponents);
		VectorBase(
			const int theDimensionCount,
			const std::vector<double, std::allocator<double>>& theLocalComponents);

	public:
		VectorBase(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalComponents);
		VectorBase(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalComponents);
		VectorBase(
			const int theDimensionCount,
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1) ;

		VectorBase(const VectorBase& rhs) = default;
		VectorBase& operator=(const VectorBase& rhs) = default;
		VectorBase(VectorBase&& rhs) = default;
		VectorBase& operator=(VectorBase&& rhs) = default;
		~VectorBase() override = default;

		bool operator==(const VectorBase& rhs) const;
		bool operator!=(const VectorBase& rhs) const;
		bool operator+=(const VectorBase& rhs) const;
		bool operator-=(const VectorBase& rhs) const;

		// Methods
		bool equalsGeometrically(ARGCOPY(VectorBase) theVector) const;

		double getLocalComponentX() const;
		double getLocalComponentY() const;
		double getLocalComponentZ() const;
		auto getLocalComponents() const -> std::array<double, 3>;
		double getGlobalComponentX() const;
		double getGlobalComponentY() const;
		double getGlobalComponentZ() const;
		auto getGlobalComponents() const -> std::array<double, 3>;
		auto getSlopes() const -> std::array<double, 3>;
		auto getAngles() const;
		double getMagnitude() const;
		auto getUnitVectorComponents() const -> std::array<double, 3>;
		void setReferenceCoordSystem(const std::shared_ptr<CoordSystem>& theCoordSystem, const bool theKeepGlobalComponentsSame);
		void setLocalComponentX(const double& theLocalComponentX);
		void setLocalComponentY(const double& theLocalComponentY);
		void setLocalComponentZ(const double& theLocalComponentZ);
		void setLocalComponents(const std::array<double, 3>& theLocalComponents);
		void setLocalComponents(const std::vector<double, std::allocator<double>>& theComponents);
		bool isParallel(ARGCOPY(VectorBase) theVector) const;
		bool isInTheSameDirection(ARGCOPY(VectorBase) theVector) const;
		bool isNormal(ARGCOPY(VectorBase) theVector) const;
		double calculateAngle(ARGCOPY(VectorBase) theVector) const;
		double dotProduct(ARGCOPY(VectorBase) theVector) const;
		auto crossProduct(ARGCOPY(VectorBase) theVector) const -> std::shared_ptr<Vector3D>;
		auto transformPoint(
			ARGCOPY(PointBase) thePoint,
			const double& theFactor) const
			-> std::shared_ptr<Point3D>;

		static auto createUnitVectorZ() -> std::shared_ptr<Vector3D>;
		static auto createUnitVectorZ(const std::shared_ptr<CoordSystem>& theCoordSystem) -> std::shared_ptr<Vector3D>;

	protected:
		void setMembers(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theCoordSystem,
			const std::array<double, 3>& theLocalComponents);
		void setMembers(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalComponents);
		template <typename It>
		void inspectLocalComponents(const It theBegin, const It theEnd) const;
		double calculateSlope(const double& theCoord0, const double& theCoord1) const;
		double calculateAngle(const double& theCoord0, const double& theCoord1) const;
		static double normalizeAngle(const double& theAngle);
		auto calculateComponentsFromAngles(const std::array<double, 3>& theAngles) const -> std::array<double, 3>;
		template<typename T>
		static T Clone(const T& arg);
	};
}

#endif
