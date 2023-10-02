/// <summary>
/// Circle defines a circle lying on a plane with a center point and radius.
/// 
/// INVARIANT:
///   Cluster of infinite number of points lying on a plane at the same distance to a point.
///   The plane is called the reference plane,
///   and the equi-distant-point is called the center point.
/// 
/// A 3D object currently.
/// Hence, does not have 2D and 3D child classes.
/// Later, will be defined as a reference object type and 2D and 3D child classes will be generated.
/// See ReferenceObject header file docstring for details.
/// 
/// See GeometryObject.hxx for the details about this library.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
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
	class Circle : public GeometryObject
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class Point3D;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Plane;

		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		std::shared_ptr<Plane> c_referencePlane = nullptr;
		std::shared_ptr<Point3D> c_centerPoint = nullptr;
		double c_radius = 0.;

		// Private default ctor used for cloning the object
		Circle() = default;

	public:
		// ctor / dtor / operators
		Circle(
			const std::shared_ptr<Plane>& theReferencePlane,
			const std::shared_ptr<Point3D>& theCenterPoint,
			const double& theRadius);
		Circle(
			Point3D thePoint0,  // Function needs a copy of the argument
			Point3D thePoint1,  // Function needs a copy of the argument
			Point3D thePoint2); // Function needs a copy of the argument

		Circle(const Circle& rhs) = default;
		Circle& operator=(const Circle& rhs) = default;
		Circle(Circle&& rhs) noexcept = default;
		Circle& operator=(Circle&& rhs) noexcept = default;
		~Circle() final = default;

		bool operator==(const Circle& rhs) const;
		bool operator!=(const Circle& rhs) const;
		bool operator+=(const Circle& rhs) const;
		bool operator-=(const Circle& rhs) const;
		
		// Methods
		bool is2D() const;
		bool is3D() const;
		bool equalsGeometrically(ARGCOPY(Circle) theCircle) const;

		auto getReferencePlane() const -> std::shared_ptr<Plane>;
		auto getCenterPoint() const -> std::shared_ptr<Point3D>;
		double getRadius() const;
		auto getCommonReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>;
		void setReferencePlane(const std::shared_ptr<Plane>& theReferencePlane);
		void setCenterPoint(const std::shared_ptr<Point3D>& theCenterPoint);
		void setRadius(const double& theRadius);
		bool isPointOn(ARGCOPY(Point3D) thePoint) const;
		bool isPointInvolved(ARGCOPY(Point3D) thePoint) const;
		int analyzePoint(ARGCOPY(Point3D) thePoint) const; // See method docstring

	private:
		void setMembers(
			const std::shared_ptr<Plane>& theReferencePlane,
			const std::shared_ptr<Point3D>& theCenterPoint,
			const double& theRadius);
	};
}

#endif
