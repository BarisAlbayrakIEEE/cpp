/// <summary>
/// Point3D defines a point wrt the reference CS.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// NO INVARIANT
/// 
/// Created for the tracebility and for the simplicity when switching between CSs.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Point3D_HeaderFile
#define _Point3D_HeaderFile

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
#ifndef _PointBase_HeaderFile
#include "PointBase.hxx"
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

	DEFINE_STANDARD_HANDLE(Point3D, PointBase)

	class Point3D : public PointBase
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
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
			return sizeof(Point3D);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		// ctor / dtor / operators
		Standard_EXPORT Point3D(const Point3D& rhs);
		Standard_EXPORT Point3D& operator=(const Point3D& rhs);
		Standard_EXPORT Point3D(Point3D&& rhs) noexcept;
		Standard_EXPORT Point3D& operator=(Point3D&& rhs) noexcept;
		Standard_EXPORT ~Point3D();

		Standard_EXPORT Point3D();
		Standard_EXPORT Point3D(const arrayS3& theLocalCoords);
		Standard_EXPORT Point3D(const vectorInput1D& theLocalCoords);
		Standard_EXPORT Point3D(ARGCOPY(CoordSystem) theReferenceCoordSystem);
		Standard_EXPORT Point3D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalCoords);
		Standard_EXPORT Point3D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords);
		Standard_EXPORT Point3D(ARGCOPY(Point2D) thePoint);

	public:
		DEFINE_STANDARD_RTTIEXT(Point3D, PointBase)
	};
}

#endif
