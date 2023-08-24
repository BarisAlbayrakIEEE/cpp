/// <summary>
/// Axis defines a directed infinite line.
/// Cluster of infinite number of points created by translating a point by a vector.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// INVARIANT:
///   The translation vector is called the direction vector,
///   and the translated point is called the passing point.
/// 
/// Although an axis is infinite, it has a direction.
/// 
/// A 3D object by default.
/// Hence, does not have 2D and 3D child classes.
/// See GeometryObject.hxx for the details.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Axis_HeaderFile
#define _Axis_HeaderFile

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

	DEFINE_STANDARD_HANDLE(Axis, GeometryObject)

	class Axis : public virtual GeometryAbstractObject, public GeometryObject
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Line;
		friend class Circle;
		friend class Plane;

	public:
		/*
		 * DEF: Modifications on OCCT
		*/
		void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(Axis);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT bool is2D() const;
		Standard_EXPORT bool is3D() const;
		Standard_EXPORT void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(Axis) theAxis) const;

	private:
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		// The axis Equation Coefficients(EC) of the form : (x - x0) / Vx = (y - y0) / Vy = (z - z0) / Vz
		// x0, y0 and z0 are coordinates of the passing point
		// Vx, Vy and Vz are the components of the direction vector
		// c_EC[i][O] = Point coord in ith direction
		// c_EC[i][1] = Vector component in ith direction
		Handle(PointBase) c_passingPoint;
		Handle(VectorBase) c_directionVector;
		arrayS32 c_EC;

	public:
		// ctor / dtor / operators
		Standard_EXPORT Axis(
			ARGCOPY(PointBase) thePassingPoint,
			ARGCOPY(VectorBase) theDirectionVector);
		Standard_EXPORT Axis(
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1);
		Standard_EXPORT Axis(
			ARGCOPY(Point2D) thePassingPoint,
			const double& theAngle);
		Standard_EXPORT Axis(
			ARGCOPY(PointBase) thePassingPoint,
			const arrayS3& theAngles);
		Standard_EXPORT Axis(const arrayS32& theEC);
		Standard_EXPORT Axis(vectorInput2D theEC);

		Standard_EXPORT Axis(const Axis& rhs);
		Standard_EXPORT Axis& operator=(const Axis& rhs);
		Standard_EXPORT Axis(Axis&& rhs) noexcept;
		Standard_EXPORT Axis& operator=(Axis&& rhs) noexcept;
		Standard_EXPORT bool operator==(const Axis& rhs);
		Standard_EXPORT bool operator!=(const Axis& rhs);
		Standard_EXPORT bool operator+=(const Axis& rhs);
		Standard_EXPORT bool operator-=(const Axis& rhs);
		Standard_EXPORT ~Axis();

	private:
		Standard_EXPORT void copyBase(const Axis& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(Axis, GeometryObject)

	public:
		Standard_EXPORT OUTVAL(PointBase) getPassingPoint() const;
		Standard_EXPORT OUTVAL(VectorBase) getDirectionVector() const;
		Standard_EXPORT arrayS32 getEC() const;
		Standard_EXPORT OUTVAL(PointBase) createPoint(const double& theFactor) const;
		Standard_EXPORT OUTVAL(Vector3D) getPassingPointAsVector() const;
		Standard_EXPORT OUTVAL(Vector3D) getPointAsVector(const double& theFactor) const;
		Standard_EXPORT OUTVAL(CoordSystem) getCommonReferenceCoordSystem() const;
		Standard_EXPORT void setPassingPoint(ARGCOPY(PointBase) thePassingPoint);
		Standard_EXPORT void setDirectionVector(ARGCOPY(VectorBase) theDirectionVector);
		Standard_EXPORT void setEC(const arrayS32& theEC);
		Standard_EXPORT void setEC(const vectorInput2D& theEC);
		Standard_EXPORT bool isParallel(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT bool isParallel(
			ARGCOPY(Axis) theAxis,
			const double& theTolerance) const;
		Standard_EXPORT bool isInTheSameDirection(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT bool isInTheSameDirection(
			ARGCOPY(Axis) theAxis,
			const double& theTolerance) const;
		Standard_EXPORT bool isNormal(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT bool isNormal(
			ARGCOPY(Axis) theAxis,
			const double& theTolerance) const;
		Standard_EXPORT bool isSkew(ARGCOPY(Axis) theAxis);
		Standard_EXPORT bool isSkew(ARGCOPY(Line) theLine);

		Standard_EXPORT bool includes(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT bool intersects(ARGCOPY(Axis) theAxis);
		Standard_EXPORT bool intersects(ARGCOPY(Line) theLine);
		Standard_EXPORT bool coincides(ARGCOPY(Axis) theAxis);
		Standard_EXPORT bool coincides(ARGCOPY(Line) theLine);
		Standard_EXPORT std::pair<INTERSECTION1, Handle(PointBase)> intersect(ARGCOPY(Axis) theAxis);
		Standard_EXPORT std::pair<INTERSECTION1, Handle(PointBase)> intersect(ARGCOPY(Line) theLine);
		Standard_EXPORT OUTVAL(PointBase) project(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT double calculateDistance(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT double calculateDistance(ARGCOPY(Axis) theAxis);
		Standard_EXPORT double calculateDistance(ARGCOPY(Line) theLine);
		Standard_EXPORT std::vector<Handle(PointBase)> findClosestPoints(ARGCOPY(Axis) theAxis);

		Standard_EXPORT double calculateCoordX_fromCoordY(const double& theCoordY) const;
		Standard_EXPORT double calculateCoordX_fromCoordZ(const double& theCoordZ) const;
		Standard_EXPORT double calculateCoordY_fromCoordX(const double& theCoordX) const;
		Standard_EXPORT double calculateCoordY_fromCoordZ(const double& theCoordZ) const;
		Standard_EXPORT double calculateCoordZ_fromCoordX(const double& theCoordX) const;
		Standard_EXPORT double calculateCoordZ_fromCoordY(const double& theCoordY) const;

	private:
		Standard_EXPORT void setMembers(const arrayS32& theEC);
		Standard_EXPORT void setMembers(const vectorInput2D& theEC);
		Standard_EXPORT void setMembers(
			ARGCOPY(PointBase) thePassingPoint,
			ARGCOPY(VectorBase) theDirectionVector);
		Standard_EXPORT void inspectEC(const arrayS32& theEC);
		Standard_EXPORT void inspectEC(const vectorInput2D& theEC);
		Standard_EXPORT int determineDimensionCountUsingEC(const arrayS32& theEC) const;
		Standard_EXPORT int determineDimensionCountUsingEC(const vectorInput2D& theEC) const;
		Standard_EXPORT void updateEC();
		Standard_EXPORT void applyEC();

		Standard_EXPORT double calculateCoord(
			int theIndexCoordO,
			int theIndexCoordl,
			const double& theCoord) const;
		Standard_EXPORT std::pair<INTERSECTION1, Handle(PointBase)> intersectBase(
			ARGCOPY(Axis) theAxis,
			const arrayS3& theCrossProductComponents) const;
		Standard_EXPORT bool equalsBase(ARGCOPY(Axis) theAxis) const;
	};
}

#endif
