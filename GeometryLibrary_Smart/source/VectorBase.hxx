/// <summary>
/// VectorBase defines a vector and its members.
/// 
/// INVARIANT:
///   A vector shall have a non-zero magnitude.
///   Hence, at least one component shall be non-zero.
/// 
/// Vector2D and Vector3D are the child classes.
/// The copy/move ctors and operators are removed to prevent slicing.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
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

	class VectorBase : public ReferenceObject
	{
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

	public:
		bool equals(ARGCOPY(VectorBase) theVector) const;
		bool equalsGeometrically(ARGCOPY(VectorBase) theVector) const;

	protected:
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		arrayS3 c_localComponents;
		double c_magnitude;

		// ctor / dtor / operators
	protected:
		VectorBase(const int theDimensionCount);
		VectorBase(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem);

	public:
		VectorBase(
			const int theDimensionCount,
			CoordSystem&theReferenceCoordSystem,
			const arrayS3& theLocalComponents);
		VectorBase(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theLocalComponents);
		VectorBase(
			const int theDimensionCount,
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1) ;

		VectorBase(const VectorBase& rhs) = delete;
		VectorBase& operator=(const VectorBase& rhs) = delete;
		VectorBase(VectorBase&& rhs) = delete;
		VectorBase& operator=(VectorBase&& rhs) = delete;
		bool operator==(const VectorBase& rhs);
		bool operator!=(const VectorBase& rhs);
		bool operator+=(const VectorBase& rhs);
		bool operator-=(const VectorBase& rhs);
		VectorBase operator+(ARGCOPY(VectorBase) theVector);
		VectorBase operator-(ARGCOPY(VectorBase) theVector);
		virtual ~VectorBase();

	private:
		void copyBase(const VectorBase& rhs);

	protected:
		void Destroy();
		static VectorBase clone(ARGCOPY(VectorBase) theVector);

	protected:
		static VectorBase DownCast(ARGCOPY(Vector2D) theVector);
		static VectorBase DownCast(ARGCOPY(Vector3D) theVector);

		// Methods
	public:
		double getLocalComponentX() const;
		double getLocalComponentY() const;
		double getLocalComponentZ() const;
		arrayS3 getLocalComponents() const;
		double getGlobalComponentX() const;
		double getGlobalComponentY() const;
		double getGlobalComponentZ() const;
		arrayS3 getGlobalComponents() const;
		arrayS3 getSlopes() const;
		arrayS3 getAngles() const;
		double getMagnitude() const;
		arrayS3 getUnitVectorComponents() const;
		void setReferenceCoordSystem(CoordSystem& theCoordSystem, const bool theKeepGlobalComponentsSame);
		void setLocalComponentX(const double& theLocalComponentX);
		void setLocalComponentY(const double& theLocalComponentY);
		void setLocalComponentZ(const double& theLocalComponentZ);
		void setLocalComponents(const arrayS3& theLocalComponents);
		void setLocalComponents(const vectorInput1D& theComponents);
		bool isParallel(ARGCOPY(VectorBase) theVector) const;
		bool isParallel(ARGCOPY(VectorBase) theVector, const double& theTolerance) const;
		bool isInTheSameDirection(ARGCOPY(VectorBase) theVector) const;
		bool isInTheSameDirection(ARGCOPY(VectorBase) theVector, const double& theTolerance) const;
		bool isNormal(ARGCOPY(VectorBase) theVector) const;
		bool isNormal(ARGCOPY(VectorBase) theVector, const double& theTolerance) const;
		double calculateAngle(ARGCOPY(VectorBase) theVector) const;
		double dotProduct(ARGCOPY(VectorBase) theVector) const;
		Vector3D crossProduct(ARGCOPY(VectorBase) theVector) const;
		Point3D transformPoint(ARGCOPY(PointBase) thePoint, const double& theFactor) const;
		VectorBase add(ARGCOPY(VectorBase) theVector) const;
		VectorBase subtruct(ARGCOPY(VectorBase) theVector) const;
		VectorBase multiply(const double& theFactor) const;

	public:
		static VectorBase createUnitVectorX(const int theDimensionCount);
		static VectorBase createUnitVectorY(const int theDimensionCount);
		static Vector3D createUnitVectorZ();
		static VectorBase createUnitVectorX(const int theDimensionCount, CoordSystem& theCoordSystem);
		static VectorBase createUnitVectorY(const int theDimensionCount, CoordSystem& theCoordSystem);
		static Vector3D createUnitVectorZ(CoordSystem& theCoordSystem);

	protected:
		void setMembers(
			const int theDimensionCount,
			CoordSystem& theCoordSystem,
			const arrayS3& theLocalComponents);
		void setMembers(
			const int theDimensionCount,
			CoordSystem& theCoordSystem,
			const vectorInput1D& theLocalComponents);
		void inspectLocalComponents(const arrayS3& theLocalComponents) const;
		double calculateMagnitude(const arrayS3& theLocalComponents);
		double calculateSlope(const double& theCoord0, const double& theCoord1) const;
		double calculateAngle(const double& theCoord0, const double& theCoord1) const;
		arrayS3 calculateComponentsFromAngles(const arrayS3& theAngles) const;

	private:
		bool equalsBase(ARGCOPY(VectorBase) theVector) const;
	};

	VectorBase operator+(ARGCOPY(VectorBase) theVector1, ARGCOPY(VectorBase) theVector2);
	VectorBase operator-(ARGCOPY(VectorBase) theVector1, ARGCOPY(VectorBase) theVector2);
}

#endif
