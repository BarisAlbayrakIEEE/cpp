/// <summary>
/// PointBase defines a point and its members.
/// 
/// NO INVARIANT
/// 
/// Point2D and Point3D are the child classes.
/// No slicing issue as the child classes do not have additional members (for tracebility and CS switcing).
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
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
	class GeometryObject;
	class ReferenceObject;
	class CoordSystem;
	class GlobalCoordSystem;
	class PointBase;
	class Point2D;
	class Point3D;
	class VectorBase;
	class Vector2D;
	class Vector3D;
	class Axis;
	class Line;
	class Circle;
	class Plane;

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

	public:
		bool equals(ARGCOPY(PointBase) thePoint) const;
		bool equalsGeometrically(ARGCOPY(PointBase) thePoint) const;

	protected:
		// Members
		arrayS3 c_localCoords;

		// ctor / dtor / operators
	protected:
		PointBase(const int theDimensionCount);
		PointBase(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem);

	public:
		PointBase(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem,
			const arrayS3& theLocalCoords);
		PointBase(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords);

		PointBase(const PointBase& rhs);
		PointBase& operator=(const PointBase& rhs);
		PointBase(PointBase&& rhs) noexcept;
		PointBase& operator=(PointBase&& rhs) noexcept;
		bool operator==(const PointBase& rhs);
		bool operator!=(const PointBase& rhs);
		bool operator+=(const PointBase& rhs);
		bool operator-=(const PointBase& rhs);
		virtual ~PointBase();

	private:
		void copyBase(const PointBase& rhs);

	protected:
		static PointBase DownCast(ARGCOPY(PointBase) thePoint);

		// Methods
	public:
		double getLocalCoordX() const;
		double getLocalCoordY() const;
		double getLocalCoordZ() const;
		arrayS3 getLocalCoords() const;
		double getGlobalCoordX() const;
		double getGlobalCoordY() const;
		double getGlobalCoordZ() const;
		arrayS3 getGlobalCoords() const;
		void setReferenceCoordSystem(
			CoordSystem& theReferenceCoordSystem,
			const bool theKeepGlobalCoordsSame);
		void setLocalCoordX(const double& theLocalCoordX);
		void setLocalCoordY(const double& theLocalCoordY);
		void setLocalCoordZ(const double& theLocalCoordZ);
		void setLocalCoords(const arrayS3& theLocalCoords);
		void setLocalCoords(const vectorInput1D& theLocalCoords);
		bool coincides(ARGCOPY(PointBase) thePoint) const;
		double calculateDistance(ARGCOPY(PointBase) thePoint) const;
		PointBase createMidPoint(ARGCOPY(PointBase) thePoint) const;
		PointBase createInterpolationPoint(
			ARGCOPY(PointBase) thePoint,
			const double& theFactor) const;

	public:
		static PointBase createPointAtOrigin(int theDimensionCount);
		static PointBase createPointAtOrigin(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem);

	protected:
		void Destroy();
		void setMembers(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem,
			const arrayS3& theLocalCoords);
		void setMembers(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords);
		static arrayS3 interpolateCoords(
			const arrayS3& theCoordsO,
			const arrayS3& theCoordsl,
			const double& theFactor);

	private:
		bool equalsBase(ARGCOPY(PointBase) thePoint) const;
	};
}

#endif
