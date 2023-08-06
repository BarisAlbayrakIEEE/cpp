/// <summary>
/// Point2D defines a point having z-coordinate equal to zero wrt the reference CS.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// INVARIANT:
///   The local coordinate in z-direction wrt the reference CS shall be zero.
/// 
/// Created for the tracebility and for the simplicity when switching to a local CS.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Point2D_HeaderFile
#define _Point2D_HeaderFile

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

	DEFINE_STANDARD_HANDLE(Point2D, PointBase)

	class Point2D : public PointBase
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class PointBase;
		friend class Point3D;
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
			return sizeof(Point2D);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		// ctor / dtor / operators
		Standard_EXPORT Point2D(const Point2D& rhs);
		Standard_EXPORT Point2D& operator=(const Point2D& rhs);
		Standard_EXPORT Point2D(Point2D&& rhs) noexcept;
		Standard_EXPORT Point2D& operator=(Point2D&& rhs) noexcept;
		Standard_EXPORT ~Point2D();

		Standard_EXPORT Point2D();
		Standard_EXPORT Point2D(const arrayS2& theLocalCoords);
		Standard_EXPORT Point2D(const arrayS3& theLocalCoords);
		Standard_EXPORT Point2D(const vectorInput1D& theLocalCoords);
		Standard_EXPORT Point2D(ARGCOPY(CoordSystem) theReferenceCoordSystem) throw (NullptrException);
		Standard_EXPORT Point2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS2& theLocalCoords) throw (NullptrException);
		Standard_EXPORT Point2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalCoords) throw (NullptrException);
		Standard_EXPORT Point2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalCoords) throw (NullptrException);

	public:
		DEFINE_STANDARD_RTTIEXT(Point2D, PointBase)

	public:
		Standard_EXPORT void setLocalCoordZ(double theLocalCoordZ) = delete;
	};
}

#endif
