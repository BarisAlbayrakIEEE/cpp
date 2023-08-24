/// <summary>
/// Vector3D defines a vector wrt the reference CS.
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
		// Members

		// ctor / dtor / operators
		Vector3D(const Vector3D& rhs);
		Vector3D& operator=(const Vector3D& rhs);
		Vector3D(Vector3D&& rhs) noexcept;
		Vector3D& operator=(Vector3D&& rhs) noexcept;
		~Vector3D();

		Vector3D(const arrayS3& theLocalComponents);
		Vector3D(const vectorInput1D& theLocalComponents);
		Vector3D(const arrayS3& theAngles, int* theNone);
		Vector3D(const vectorInput1D& theAngles, int* theNone);
		Vector3D(ARGCOPY(PointBase) thePoint);
		Vector3D(
			ARGCOPY(PointBase) thePoint0,
			ARGCOPY(PointBase) thePoint1);
		Vector3D(
			CoordSystem& theReferenceCoordSystem,
			const arrayS3& theLocalComponents);
		Vector3D(
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theLocalComponents);
		Vector3D(
			CoordSystem& theReferenceCoordSystem,
			const arrayS3& theAngles,
			int* theNone);
		Vector3D(
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theAngles,
			int* theNone);
		Vector3D(ARGCOPY(Vector2D) theVector);

		// Methods
	public:
		Vector3D getUnitVector() const;
		void setLocalComponentsUsingAngles(const arrayS3& theAngles);
		void setLocalComponentsUsingAngles(const vectorInput1D& theAngles);
	};
}

#endif
