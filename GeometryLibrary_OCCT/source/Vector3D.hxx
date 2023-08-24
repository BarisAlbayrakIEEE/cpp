/// <summary>
/// Vector3D defines a vector wrt the reference CS.
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
/// Created for the tracebility and for the simplicity when switching between CSs.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Vector3D_HeaderFile
#define _Vector3D_HeaderFile

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

	DEFINE_STANDARD_HANDLE(Vector3D, VectorBase)

	class Vector3D : public VectorBase
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
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
			return sizeof(Vector3D);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		// Members

		// ctor / dtor / operators
		Standard_EXPORT Vector3D(const Vector3D& rhs);
		Standard_EXPORT Vector3D& operator=(const Vector3D& rhs);
		Standard_EXPORT Vector3D(Vector3D&& rhs) noexcept;
		Standard_EXPORT Vector3D& operator=(Vector3D&& rhs) noexcept;
		Standard_EXPORT ~Vector3D();

		Standard_EXPORT Vector3D(const arrayS3& theLocalComponents);
		Standard_EXPORT Vector3D(const vectorInput1D& theLocalComponents);
		Standard_EXPORT Vector3D(const arrayS3& theAngles, int* theNone);
		Standard_EXPORT Vector3D(const vectorInput1D& theAngles, int* theNone);
		Standard_EXPORT Vector3D(ARGCOPY(PointBase) thePoint);
		Standard_EXPORT Vector3D(
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1);
		Standard_EXPORT Vector3D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theLocalComponents);
		Standard_EXPORT Vector3D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theLocalComponents);
		Standard_EXPORT Vector3D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const arrayS3& theAngles,
			int* theNone);
		Standard_EXPORT Vector3D(
			ARGCOPY(CoordSystem) theReferenceCoordSystem,
			const vectorInput1D& theAngles,
			int* theNone);
		Standard_EXPORT Vector3D(ARGCOPY(Vector2D) theVector);

		// Methods
	public:
		Standard_EXPORT OUTVAL(Vector3D) getUnitVector() const;
		Standard_EXPORT void setLocalComponentsUsingAngles(const arrayS3& theAngles);
		Standard_EXPORT void setLocalComponentsUsingAngles(const vectorInput1D& theAngles);

	public:
		DEFINE_STANDARD_RTTIEXT(Vector3D, VectorBase)

	public:
	};
}

#endif
