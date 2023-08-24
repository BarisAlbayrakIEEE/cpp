/// <summary>
/// Vector2D defines a vector having z-component equal to zero wrt the reference CS.
/// 
/// INVARIANT:
///   1. A vector shall have a non-zero magnitude. Hence, at least one component shall be non-zero.
///   2. The local component in z-direction wrt the reference CS shall be zero.
/// 
/// Created for the tracebility and for the simplicity when switching to a local CS.
/// 
/// See GeometryObject.hxx for details
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Vector2D_HeaderFile
#define _Vector2D_HeaderFile

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
		// Members

		// ctor / dtor / operators
		Vector2D(const Vector2D& rhs);
		Vector2D& operator=(const Vector2D& rhs);
		Vector2D(Vector2D&& rhs) noexcept;
		Vector2D& operator=(Vector2D&& rhs) noexcept;
		~Vector2D();

		Vector2D(const arrayS2& theLocalComponents);
		Vector2D(const arrayS3& theLocalComponents);
		Vector2D(const vectorInput1D& theLocalComponents);
		Vector2D(const double& theAngle);
		Vector2D(ARGCOPY(Point2D) thePoint);
		Vector2D(
			ARGCOPY(Point2D) thePoint0,
			ARGCOPY(Point2D) thePoint1);
		Vector2D(
			CoordSystem& theReferenceCoordSystem,
			const arrayS2& theLocalComponents);
		Vector2D(
			CoordSystem& theReferenceCoordSystem,
			const arrayS3& theLocalComponents);
		Vector2D(
			CoordSystem& theReferenceCoordSystem,
			const vectorInput1D& theLocalComponents);
		Vector2D(
			CoordSystem& theReferenceCoordSystem,
			const double& theAngle);

		// Methods
	public:
		Vector2D getUnitVector() const;
		double getSlope();
		double getAngle();
		Vector2D createNormalVector();
		void setLocalComponentZ(const double& theLocalComponentZ) = delete;
	};
}

#endif
