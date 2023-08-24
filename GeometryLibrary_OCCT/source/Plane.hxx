/// <summary>
/// Plane defines a plane with a passing point and normal vector.
/// See GeometryObject.hxx for details
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// NO INVARIANT
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Plane_HeaderFile
#define _Plane_HeaderFile

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

	DEFINE_STANDARD_HANDLE(Plane, GeometryObject)

	class Plane : public virtual GeometryAbstractObject, public GeometryObject
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;

	public:
		/*
		 * DEF: Modifications on OCCT
		*/
		void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(Plane);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT bool is2D() const;
		Standard_EXPORT bool is3D() const;
		Standard_EXPORT void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(Plane) thePlane) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(Plane) thePlane) const;

	private:
		// Members
		Handle(PointBase) c_passingPoint;
		Handle(VectorBase) c_normalVector;
		arrayS4 c_EC; // C0x + C1y + C2z + C3 = O

	public:
		// ctor / dtor / operators
		Standard_EXPORT Plane(const arrayS4& theEC);
		Standard_EXPORT Plane(const vectorInput1D& theEC);
		Standard_EXPORT Plane(
			ARGCOPY(PointBase) thePoint,
			ARGCOPY(VectorBase) theNormalVector);
		Standard_EXPORT Plane(
			ARGCOPY(PointBase) thePoint,
			ARGCOPY(VectorBase) theInPlaneVectorO,
			ARGCOPY(VectorBase) theInPlaneVector1);
		Standard_EXPORT Plane(
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1,
			ARGCOPY(PointBase) thePoint2);

		Standard_EXPORT Plane(const Plane& rhs);
		Standard_EXPORT Plane& operator=(const Plane& rhs);
		Standard_EXPORT Plane(Plane&& rhs) noexcept;
		Standard_EXPORT Plane& operator=(Plane&& rhs) noexcept;
		Standard_EXPORT bool operator==(const Plane& rhs);
		Standard_EXPORT bool operator!=(const Plane& rhs);
		Standard_EXPORT bool operator+=(const Plane& rhs);
		Standard_EXPORT bool operator-=(const Plane& rhs);
		Standard_EXPORT ~Plane();

	private:
		Standard_EXPORT void copyBase(const Plane& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(Plane, GeometryObject)

			// Methods
	public:
		Standard_EXPORT OUTVAL(PointBase) getPassingPoint() const;
		Standard_EXPORT OUTVAL(VectorBase) getNormalVector() const;
		Standard_EXPORT arrayS4 getEC() const;
		Standard_EXPORT OUTVAL(CoordSystem) getCommonReferenceCoordSystem() const;
		Standard_EXPORT void setMembers(const arrayS4& theEC);
		Standard_EXPORT void setMembers(const vectorInput1D& theEC);
		Standard_EXPORT void setMembers(
			ARGCOPY(PointBase) thePassingPoint,
			ARGCOPY(VectorBase) theNormalVector);
		Standard_EXPORT void setPassingPoint(ARGCOPY(PointBase) thePassingPoint);
		Standard_EXPORT void setNormalVector(ARGCOPY(VectorBase) theNormalVector);
		Standard_EXPORT void setEC(const arrayS4& theEC);
		Standard_EXPORT void setEC(const vectorInput1D& theEC);

		Standard_EXPORT bool intersects(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT bool intersects(ARGCOPY(Line) theLine) const;
		Standard_EXPORT bool includes(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT bool includes(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT bool includes(ARGCOPY(Line) theLine) const;
		Standard_EXPORT std::pair<bool, double> calculateCoordX(const double& theCoordY, const double& theCoordZ) const;
		Standard_EXPORT std::pair<bool, double> calculateCoordY(const double& theCoordZ, const double& theCoordX) const;
		Standard_EXPORT std::pair<bool, double> calculateCoordZ(const double& theCoordX, const double& theCoordY) const;
		Standard_EXPORT std::pair<bool, std::array<double, 2>> calculateCoordXY(const double& theCoordZ) const;
		Standard_EXPORT std::pair<bool, std::array<double, 2>> calculateCoordYZ(const double& theCoordX) const;
		Standard_EXPORT std::pair<bool, std::array<double, 2>> calculateCoordZX(const double& theCoordY) const;
		Standard_EXPORT double calculateDistance(ARGCOPY(Plane) thePlane) const;
		Standard_EXPORT double calculateDistance(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT double calculateDistance(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT double calculateDistance(ARGCOPY(Line) theLine) const;
		Standard_EXPORT std::pair<INTERSECTION2, Handle(Axis)> intersect(ARGCOPY(Plane) thePlane) const;
		Standard_EXPORT std::pair<INTERSECTION2, Handle(Point3D)> intersect(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT std::pair<INTERSECTION2, Handle(Point3D)> intersect(ARGCOPY(Line) theLine) const;
		Standard_EXPORT OUTVAL(Point3D) project(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT OUTVAL(Vector3D) project(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT OUTVAL(Axis) project(ARGCOPY(Axis) theAxis) const;
		Standard_EXPORT OUTVAL(Line) project(ARGCOPY(Line) theLine) const;
		Standard_EXPORT OUTVAL(Point3D) createPoint(const double& theFactor) const;

	private:
		Standard_EXPORT bool equalsBase(ARGCOPY(Plane) thePlane) const;
		Standard_EXPORT void inspectEC(const arrayS4& theEC);
		Standard_EXPORT void inspectEC(const vectorInput1D& theEC);
		Standard_EXPORT void updateEC();
		Standard_EXPORT void applyEC();
	};
}

#endif
