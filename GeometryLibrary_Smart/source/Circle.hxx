/// <summary>
/// Circle defines a circle lying on a plane with a center point and radius.
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

	class Circle : public GeometryObject
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
		bool is2D() const;
		bool is3D() const;

	public:
		bool equals(ARGCOPY(Circle) theCircle) const;
		bool equalsGeometrically(ARGCOPY(Circle) theCircle) const;

	private:
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		std::shared_ptr<Plane> c_referencePlane;
		std::shared_ptr<PointBase> c_centerPoint;
		double c_radius;

	public:
		// ctor / dtor / operators
		Circle(
			Plane& theReferencePlane,
			PointBase& theCenterPoint,
			const double& theRadius);
		Circle(
			PointBase& thePoint0,
			PointBase& thePoint1,
			PointBase& thePoint2);

		Circle(const Circle& rhs);
		Circle& operator=(const Circle& rhs);
		Circle(Circle&& rhs) noexcept;
		Circle& operator=(Circle&& rhs) noexcept;
		bool operator==(const Circle& rhs);
		bool operator!=(const Circle& rhs);
		bool operator+=(const Circle& rhs);
		bool operator-=(const Circle& rhs);
		~Circle();

	private:
		void copyBase(const Circle& rhs);
		void Destroy();
		
		// Methods
	public:
		Plane& getReferencePlane() const;
		PointBase& getCenterPoint() const;
		double getRadius() const;
		CoordSystem* getCommonReferenceCoordSystem() const;
		void setReferencePlane(Plane& theReferencePlane);
		void setCenterPoint(PointBase& theCenterPoint);
		void setRadius(const double& theRadius);
		bool isPointOn(ARGCOPY(PointBase) thePoint) const;
		bool isPointInvolved(ARGCOPY(PointBase) thePoint) const;
		int analyzePoint(ARGCOPY(PointBase) thePoint) const; // See method docstring

	private:
		void setMembers(
			Plane& theReferencePlane,
			PointBase& theCenterPoint,
			const double& theRadius);

		bool equalsBase(ARGCOPY(Circle) theCircle) const;
	};
}

#endif
