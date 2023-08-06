/// <summary>
/// Vector2D defines a vector having z-component equal to zero wrt the reference CS.
/// 
/// INVARIANT:
///   1. A vector shall have a non-zero magnitude. Hence, at least one component shall be non-zero.
///   2. The local component in z-direction wrt the reference CS shall be zero.
/// 
/// Created for the tracebility and for the simplicity when switching to a local CS.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// NO INVARIANT
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Vector2D_HeaderFile
#define _Vector2D_HeaderFile

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
#ifndef _VectorBase_HeaderFile
#include "VectorBase.hxx"
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

	DEFINE_STANDARD_HANDLE(Vector2D, VectorBase)

	class Vector2D : public VectorBase
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
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
			return sizeof(Vector2D);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		// Members

		// ctor / dtor / operators
		Standard_EXPORT Vector2D(const Vector2D& rhs);
		Standard_EXPORT Vector2D& operator=(const Vector2D& rhs);
		Standard_EXPORT Vector2D(Vector2D&& rhs) noexcept;
		Standard_EXPORT Vector2D& operator=(Vector2D&& rhs) noexcept;
		Standard_EXPORT ~Vector2D();

		Standard_EXPORT Vector2D(const arrayS2& theLocalComponents) throw (ZeroVectorException);
		Standard_EXPORT Vector2D(const arrayS3& theLocalComponents) throw (ZeroVectorException);
		Standard_EXPORT Vector2D(const vectorInput1D& theLocalComponents) throw (ZeroVectorException);
		Standard_EXPORT Vector2D(const double& theAngle);
		Standard_EXPORT Vector2D(ARGCOPY(Point2D) thePoint) throw (NullptrException, ZeroVectorException);
		Standard_EXPORT Vector2D(
			ARGCOPY(Point2D) thePoint0,
			ARGCOPY(Point2D) thePoint1) throw (NullptrException, CoordSystemMismatchException, ZeroVectorException);
		Standard_EXPORT Vector2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS2& theLocalComponents) throw (NullptrException, ZeroVectorException);
		Standard_EXPORT Vector2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalComponents) throw (NullptrException, ZeroVectorException);
		Standard_EXPORT Vector2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalComponents) throw (NullptrException, ZeroVectorException);
		Standard_EXPORT Vector2D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const double& theAngle) throw (NullptrException);

	public:
		DEFINE_STANDARD_RTTIEXT(Vector2D, VectorBase)

		// Methods
	public:
		Standard_EXPORT OUTVAL(Vector2D) getUnitVector() const;
		Standard_EXPORT double getSlope();
		Standard_EXPORT double getAngle();
		Standard_EXPORT OUTVAL(Vector2D) createNormalVector();
		Standard_EXPORT void setLocalComponentZ(const double& theLocalComponentZ) = delete;
	};
}

#endif
