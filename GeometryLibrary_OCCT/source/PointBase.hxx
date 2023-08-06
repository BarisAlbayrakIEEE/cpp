/// <summary>
/// PointBase defines a point and its members.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
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

	DEFINE_STANDARD_HANDLE(PointBase, ReferenceObject)

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
		/*
		 * DEF: Modifications on OCCT
		*/
		void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(PointBase);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(PointBase) thePoint) const;

	protected:
		// Members
		arrayS3 c_localCoords;

		// ctor / dtor / operators
	protected:
		Standard_EXPORT PointBase(const int theDimensionCount);
		Standard_EXPORT PointBase(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem) throw (NullptrException);

	public:
		Standard_EXPORT PointBase(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalCoords);
		Standard_EXPORT PointBase(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords);

		Standard_EXPORT PointBase(const PointBase& rhs);
		Standard_EXPORT PointBase& operator=(const PointBase& rhs);
		Standard_EXPORT PointBase(PointBase&& rhs) noexcept;
		Standard_EXPORT PointBase& operator=(PointBase&& rhs) noexcept;
		Standard_EXPORT bool operator==(const PointBase& rhs);
		Standard_EXPORT bool operator!=(const PointBase& rhs);
		Standard_EXPORT bool operator+=(const PointBase& rhs);
		Standard_EXPORT bool operator-=(const PointBase& rhs);
		Standard_EXPORT virtual ~PointBase();

	private:
		Standard_EXPORT void copyBase(const PointBase& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(PointBase, ReferenceObject)

		// Methods
	public:
		Standard_EXPORT double getLocalCoordX() const;
		Standard_EXPORT double getLocalCoordY() const;
		Standard_EXPORT double getLocalCoordZ() const;
		Standard_EXPORT arrayS3 getLocalCoords() const;
		Standard_EXPORT double getGlobalCoordX() const;
		Standard_EXPORT double getGlobalCoordY() const;
		Standard_EXPORT double getGlobalCoordZ() const;
		Standard_EXPORT arrayS3 getGlobalCoords() const;
		Standard_EXPORT void setReferenceCoordSystem(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const bool theKeepGlobalCoordsSarne);
		Standard_EXPORT void setLocalCoordX(const double& theLocalCoordX);
		Standard_EXPORT void setLocalCoordY(const double& theLocalCoordY);
		Standard_EXPORT void setLocalCoordZ(const double& theLocalCoordZ);
		Standard_EXPORT void setLocalCoords(const arrayS3& theLocalCoords);
		Standard_EXPORT void setLocalCoords(const vectorInput1D& theLocalCoords);
		Standard_EXPORT bool coincides(ARGCOPY(PointBase) thePoint) const throw (NullptrException);
		Standard_EXPORT double calculateDistance(ARGCOPY(PointBase) thePoint) const throw (NullptrException);
		Standard_EXPORT OUTVAL(PointBase) createMidPoint(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT OUTVAL(PointBase) createInterpolationPoint(
			ARGCOPY(PointBase) thePoint,
			const double& theFactor) const throw (NullptrException);

	public:
		Standard_EXPORT static OUTVAL(PointBase) createPointAtOrigin(int theDimensionCount);
		Standard_EXPORT static OUTVAL(PointBase) createPointAtOrigin(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem);

	protected:
		Standard_EXPORT void setMembers(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalCoords);
		Standard_EXPORT void setMembers(
			const int theDimensionCount,
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords);
		Standard_EXPORT static arrayS3 interpolateCoords(
			const arrayS3& theCoordsO,
			const arrayS3& theCoordsl,
			const double& theFactor);

	private:
		Standard_EXPORT bool equalsBase(ARGCOPY(PointBase) thePoint) const;
	};
}

#endif
