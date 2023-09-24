/// <summary>
/// PointBase defines a point and its members.
/// 
/// NO INVARIANT
/// 
/// Point2D and Point3D are the child classes.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4250)
#pragma warning(disable : 4290)

#ifndef _PointBase_HeaderFile
#define _PointBase_HeaderFile

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

namespace GeometryNamespace {
	class PointBase : public ReferenceObject
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

		using ReferenceObject::is2D; // To disable 4250 warning
		using ReferenceObject::is3D; // To disable 4250 warning

		// Members
		std::array<double, 3> c_localCoords = { {} };

	protected:
		// ctor / dtor / operators
		PointBase() = default;
		explicit PointBase(const int theDimensionCount);
		PointBase(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
		PointBase(
			const int theDimensionCount,
			const std::array<double, 3>& theLocalCoords);
		PointBase(
			const int theDimensionCount,
			const std::vector<double, std::allocator<double>>& theLocalCoords);

	public:
		PointBase(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalCoords);
		PointBase(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalCoords);

		PointBase(const PointBase& rhs) = default;
		PointBase& operator=(const PointBase& rhs) = default;
		PointBase(PointBase&& rhs) noexcept = default;
		PointBase& operator=(PointBase&& rhs) noexcept = default;
		~PointBase() override = default;

		bool operator==(const PointBase& rhs) const;
		bool operator!=(const PointBase& rhs) const;
		bool operator+=(const PointBase& rhs) const;
		bool operator-=(const PointBase& rhs) const;

		// Methods
		bool equals(ARGCOPY(PointBase) thePoint) const;
		bool equalsGeometrically(ARGCOPY(PointBase) thePoint) const;

		double getLocalCoordX() const;
		double getLocalCoordY() const;
		double getLocalCoordZ() const;
		auto getLocalCoords() const -> std::array<double, 3>;
		double getGlobalCoordX() const;
		double getGlobalCoordY() const;
		double getGlobalCoordZ() const;
		auto getGlobalCoords() const -> std::array<double, 3>;
		void setReferenceCoordSystem(
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const bool theKeepGlobalCoordsSame);
		void setLocalCoordX(const double& theLocalCoordX);
		void setLocalCoordY(const double& theLocalCoordY);
		void setLocalCoordZ(const double& theLocalCoordZ);
		void setLocalCoords(const std::array<double, 3>& theLocalCoords);
		void setLocalCoords(const std::vector<double, std::allocator<double>>& theLocalCoords);
		bool coincides(ARGCOPY(PointBase) thePoint) const;
		double calculateDistance(ARGCOPY(PointBase) thePoint) const;
		auto createMidPoint(ARGCOPY(PointBase) thePoint) const -> std::shared_ptr<PointBase>;
		auto createInterpolationPoint(
			ARGCOPY(PointBase) thePoint,
			const double& theFactor) const
			-> std::shared_ptr<PointBase>;

		static auto createPointAtOrigin(int theDimensionCount) -> std::shared_ptr<PointBase>;
		static auto createPointAtOrigin(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
			-> std::shared_ptr<PointBase>;

	protected:
		void setMembers(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::array<double, 3>& theLocalCoords);
		void setMembers(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
			const std::vector<double, std::allocator<double>>& theLocalCoords);
		static std::array<double, 3> interpolateCoords(
			const std::array<double, 3>& theCoordsO,
			const std::array<double, 3>& theCoordsl,
			const double& theFactor);
		template<typename T>
		static T Clone(const T& arg);
	};
}

#endif
