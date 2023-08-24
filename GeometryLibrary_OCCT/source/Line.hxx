/// <summary>
/// Line defines a line piece between two points.
/// Cluster of the points between two end points.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// INVARIANT:
///   The end points shall be geometrically unequal (i.e. non-coincident).
/// 
/// The points can be defined wrt the different refereence CSs.
/// 
/// A 3D object by default.
/// Hence, does not have 2D and 3D child classes.
/// See GeometryObject.hxx for the details.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Line_HeaderFile
#define _Line_HeaderFile

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

	DEFINE_STANDARD_HANDLE(Line, GeometryObject)

	class Line : public virtual GeometryAbstractObject, public GeometryObject
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
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
			return sizeof(Line);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT bool is2D() const;
		Standard_EXPORT bool is3D() const;
		Standard_EXPORT void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(Line) theLine) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(Line) theLine) const;

	private:
		// Members
		Handle(Axis) c_axis;
		Handle(PointBase) c_endPoint0;
		Handle(PointBase) c_endPoint1;
		double c_length;

	public:
		// ctor / dtor / operators
		Standard_EXPORT Line(
			ARGCOPY(PointBase) theEndPoint0,
			ARGCOPY(PointBase) theEndPoint1);

		Standard_EXPORT Line(const Line& rhs);
		Standard_EXPORT Line& operator=(const Line& rhs);
		Standard_EXPORT Line(Line&& rhs) noexcept;
		Standard_EXPORT Line& operator=(Line&& rhs) noexcept;
		Standard_EXPORT bool operator==(const Line& rhs);
		Standard_EXPORT bool operator!=(const Line& rhs);
		Standard_EXPORT bool operator+=(const Line& rhs);
		Standard_EXPORT bool operator-=(const Line& rhs);
		Standard_EXPORT ~Line();

	private:
		Standard_EXPORT void copyBase(const Line& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(Line, GeometryObject)

		// Methods:
	public:
		Standard_EXPORT bool includes(ARGCOPY(PointBase) thePoint);
		Standard_EXPORT bool intersects(ARGCOPY(Axis) theAxis);
		Standard_EXPORT bool intersects(ARGCOPY(Line) theLine);
		Standard_EXPORT bool coincides(ARGCOPY(Axis) theAxis);
		Standard_EXPORT bool coincides(ARGCOPY(Line) theLine);
		Standard_EXPORT bool isSkew(ARGCOPY(Axis) theAxis);
		Standard_EXPORT bool isSkew(ARGCOPY(Line) theLine);
		Standard_EXPORT std::pair<INTERSECTION1, Handle(PointBase)> intersect(ARGCOPY(Axis) theAxis);
		Standard_EXPORT std::pair<INTERSECTION1, Handle(PointBase)> intersect(ARGCOPY(Line) theLine);
		Standard_EXPORT OUTVAL(PointBase) project(ARGCOPY(PointBase) thePoint);
		Standard_EXPORT double calculateDistance(ARGCOPY(PointBase) thePoint);
		Standard_EXPORT double calculateDistance(ARGCOPY(Axis) theAxis);
		Standard_EXPORT double calculateDistance(ARGCOPY(Line) theLine);

		Standard_EXPORT OUTVAL(Axis) getAxis() const;
		Standard_EXPORT OUTVAL(VectorBase) getDirectionVector() const;
		Standard_EXPORT std::vector<Handle(PointBase)> getEndPoints() const;
		Standard_EXPORT OUTVAL(PointBase) getEndPoint0() const;
		Standard_EXPORT OUTVAL(PointBase) getEndPoint1() const;
		Standard_EXPORT double getLength() const;
		Standard_EXPORT OUTVAL(CoordSystem) getCommonReferenceCoordSystem() const;
		Standard_EXPORT void setAxis(ARGCOPY(Axis) theAxis);
		Standard_EXPORT void setEndPoint0(ARGCOPY(PointBase) theEndPoint0);
		Standard_EXPORT void setEndPoint1(ARGCOPY(PointBase) theEndPoint1);
		Standard_EXPORT OUTVAL(PointBase) createMidpoint() const;

	private:
		Standard_EXPORT bool equalsBase(ARGCOPY(Line) theLine) const;
	};
}

#endif
