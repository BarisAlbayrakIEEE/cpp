/// <summary>
/// Circle defines a circle lying on a plane with a center point and radius.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// INVARIANT:
///   Cluster of infinite number of points lying on a plane at the same distance to a point.
///   The plane is called the reference plane,
///   and the equi-distant-point is called the center point.
/// 
/// A 3D object by default.
/// Hence, does not have 2D and 3D child classes.
/// See GeometryObject.hxx for the details.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Circle_HeaderFile
#define _Circle_HeaderFile

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
#ifndef _GeometryAbstractObject_HeaderFile
#include "GeometryAbstractObject.hxx"
#endif
#ifndef _GeometryObject_HeaderFile
#include "GeometryObject.hxx"
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

	DEFINE_STANDARD_HANDLE(Circle, GeometryObject)

	class Circle : public virtual GeometryAbstractObject, public GeometryObject
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Plane;

	public:
		/*
		 * DEF: Modifications on OCCT
		*/
		void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(Circle);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT bool is2D() const;
		Standard_EXPORT bool is3D() const;
		Standard_EXPORT void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(Circle) theCircle) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(Circle) theCircle) const;

	private:
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		Handle(Plane) c_referencePlane;
		Handle(PointBase) c_centerPoint;
		double c_radius;

	public:
		// ctor / dtor / operators
		Standard_EXPORT Circle(
			ARGCOPY(Plane) theReferencePlane,
			ARGCOPY(PointBase) theCenterPoint,
			const double& theRadius) throw (
				NullptrException,
				ZeroDimensionException,
				GeometricalMismatchException);
		Standard_EXPORT Circle(
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1,
			ARGCOPY(PointBase) thePoint2) throw (
				NullptrException,
				CoincidenceException,
				ColinearPointsException,
				UncaughtException);

		Standard_EXPORT Circle(const Circle& rhs);
		Standard_EXPORT Circle& operator=(const Circle& rhs);
		Standard_EXPORT Circle(Circle&& rhs) noexcept;
		Standard_EXPORT Circle& operator=(Circle&& rhs) noexcept;
		Standard_EXPORT bool operator==(const Circle& rhs);
		Standard_EXPORT bool operator!=(const Circle& rhs);
		Standard_EXPORT bool operator+=(const Circle& rhs);
		Standard_EXPORT bool operator-=(const Circle& rhs);
		Standard_EXPORT ~Circle();

	private:
		void copyBase(const Circle& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(Circle, GeometryObject)
		
		// Methods
	public:
		Standard_EXPORT OUTVAL(Plane) getReferencePlane() const;
		Standard_EXPORT OUTVAL(PointBase) getCenterPoint() const;
		Standard_EXPORT double getRadius() const;
		Standard_EXPORT OUTVAL(CoordSystem) getCommonReferenceCoordSystem() const;
		Standard_EXPORT void setReferencePlane(ARGCOPY(Plane) theReferencePlane) throw (
			NullptrException,
			GeometricalMismatchException);
		Standard_EXPORT void setCenterPoint(ARGCOPY(PointBase) theCenterPoint) throw (
			NullptrException,
			GeometricalMismatchException);
		Standard_EXPORT void setRadius(const double& theRadius) throw (ZeroDimensionException);
		Standard_EXPORT bool isPointOn(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT bool isPointInvolved(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT int analyzePoint(ARGCOPY(PointBase) thePoint) const throw (NullptrException); // See method docstring

	private:
		Standard_EXPORT void setMembers(
			ARGCOPY(Plane) theReferencePlane,
			ARGCOPY(PointBase) theCenterPoint,
			const double& theRadius) throw (
				NullptrException,
				ZeroDimensionException);

		Standard_EXPORT bool equalsBase(ARGCOPY(Circle) theCircle) const;
	};
}

#endif
