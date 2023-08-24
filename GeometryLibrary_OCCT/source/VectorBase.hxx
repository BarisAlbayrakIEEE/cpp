/// <summary>
/// VectorBase defines a vector and its members.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
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

#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Handle_HeaderFile
#include <Standard_Handle.hxx>
#endif
#ifndef _Standard_Type_HeaderFile
#include <Standard_Type.hxx>
#endif
#ifndef _Standard_Size_HeaderFile
#include <Standard_Size.hxx>
#endif
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

	DEFINE_STANDARD_HANDLE(VectorBase, ReferenceObject)

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
		/*
		 * DEF: Modifications on OCCT
		*/
		void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(VectorBase);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(VectorBase) theVector) const;

	protected:
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		arrayS3 c_localComponents;
		double c_magnitude;

		// ctor / dtor / operators
	protected:
		Standard_EXPORT VectorBase(const int theDimensionCount);
		Standard_EXPORT VectorBase(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem);

	public:
		Standard_EXPORT VectorBase(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalComponents);
		Standard_EXPORT VectorBase(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalComponents);
		Standard_EXPORT VectorBase(
			const int theDimensionCount,
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1) ;

		Standard_EXPORT VectorBase(const VectorBase& rhs) = delete;
		Standard_EXPORT VectorBase& operator=(const VectorBase& rhs) = delete;
		Standard_EXPORT VectorBase(VectorBase&& rhs) = delete;
		Standard_EXPORT VectorBase& operator=(VectorBase&& rhs) = delete;
		Standard_EXPORT bool operator==(const VectorBase& rhs);
		Standard_EXPORT bool operator!=(const VectorBase& rhs);
		Standard_EXPORT bool operator+=(const VectorBase& rhs);
		Standard_EXPORT bool operator-=(const VectorBase& rhs);
		Standard_EXPORT OUTVAL(VectorBase) operator+(ARGCOPY(VectorBase) theVector);
		Standard_EXPORT OUTVAL(VectorBase) operator-(ARGCOPY(VectorBase) theVector);
		Standard_EXPORT virtual ~VectorBase();

	private:
		Standard_EXPORT void copyBase(const VectorBase& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(VectorBase, ReferenceObject)

		// Methods
	public:
		Standard_EXPORT double getLocalComponentX() const;
		Standard_EXPORT double getLocalComponentY() const;
		Standard_EXPORT double getLocalComponentZ() const;
		Standard_EXPORT arrayS3 getLocalComponents() const;
		Standard_EXPORT double getGlobalComponentX() const;
		Standard_EXPORT double getGlobalComponentY() const;
		Standard_EXPORT double getGlobalComponentZ() const;
		Standard_EXPORT arrayS3 getGlobalComponents() const;
		Standard_EXPORT arrayS3 getSlopes() const;
		Standard_EXPORT arrayS3 getAngles() const;
		Standard_EXPORT double getMagnitude() const;
		Standard_EXPORT arrayS3 getUnitVectorComponents() const;
		Standard_EXPORT void setReferenceCoordSystem(ARGCOPY(CoordSystem) theCoordSystem, const bool theKeepGlobalComponentsSame);
		Standard_EXPORT void setLocalComponentX(const double& theLocalComponentX);
		Standard_EXPORT void setLocalComponentY(const double& theLocalComponentY);
		Standard_EXPORT void setLocalComponentZ(const double& theLocalComponentZ);
		Standard_EXPORT void setLocalComponents(const arrayS3& theLocalComponents);
		Standard_EXPORT void setLocalComponents(const vectorInput1D& theComponents);
		Standard_EXPORT bool isParallel(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT bool isParallel(ARGCOPY(VectorBase) theVector, const double& theTolerance) const;
		Standard_EXPORT bool isInTheSameDirection(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT bool isInTheSameDirection(ARGCOPY(VectorBase) theVector, const double& theTolerance) const;
		Standard_EXPORT bool isNormal(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT bool isNormal(ARGCOPY(VectorBase) theVector, const double& theTolerance) const;
		Standard_EXPORT double calculateAngle(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT double dotProduct(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT OUTVAL(Vector3D) crossProduct(ARGCOPY(VectorBase) theVector);
		Standard_EXPORT OUTVAL(Point3D) transformPoint(ARGCOPY(PointBase) thePoint, const double& theFactor) const;
		Standard_EXPORT OUTVAL(VectorBase) add(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT OUTVAL(VectorBase) subtruct(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT OUTVAL(VectorBase) multiply(const double& theFactor) const;

	public:
		Standard_EXPORT static OUTVAL(VectorBase) createUnitVectorX(const int theDimensionCount);
		Standard_EXPORT static OUTVAL(VectorBase) createUnitVectorY(const int theDimensionCount);
		Standard_EXPORT static OUTVAL(Vector3D) createUnitVectorZ();
		Standard_EXPORT static OUTVAL(VectorBase) createUnitVectorX(const int theDimensionCount, ARGCOPY(CoordSystem) theCoordSystem);
		Standard_EXPORT static OUTVAL(VectorBase) createUnitVectorY(const int theDimensionCount, ARGCOPY(CoordSystem) theCoordSystem);
		Standard_EXPORT static OUTVAL(Vector3D) createUnitVectorZ(ARGCOPY(CoordSystem) theCoordSystem);

	protected:
		Standard_EXPORT void setMembers(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theCoordSystem,
			const arrayS3& theLocalComponents);
		Standard_EXPORT void setMembers(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theCoordSystem,
			const vectorInput1D& theLocalComponents);
		Standard_EXPORT void inspectLocalComponents(const arrayS3& theLocalComponents) const;
		Standard_EXPORT double calculateMagnitude(const arrayS3& theLocalComponents);
		Standard_EXPORT double calculateSlope(const double& theCoord0, const double& theCoord1) const;
		Standard_EXPORT double calculateAngle(const double& theCoord0, const double& theCoord1) const;
		Standard_EXPORT arrayS3 calculateComponentsFromAngles(const arrayS3& theAngles) const;

	private:
		Standard_EXPORT bool equalsBase(ARGCOPY(VectorBase) theVector) const;
	};

	Standard_EXPORT OUTVAL(VectorBase) operator+(ARGCOPY(VectorBase) theVector1, ARGCOPY(VectorBase) theVector2);
	Standard_EXPORT OUTVAL(VectorBase) operator-(ARGCOPY(VectorBase) theVector1, ARGCOPY(VectorBase) theVector2);
}

#endif
