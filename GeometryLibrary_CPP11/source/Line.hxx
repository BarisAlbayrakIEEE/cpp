/// <summary>
/// Line defines a line piece between two points.
/// Cluster of the points between two end points.
/// 
/// INVARIANT:
///   The end points shall be geometrically unequal (i.e. non-coincident).
/// 
/// The points can be defined wrt different reference CSs.
/// 
/// A 3D object currently.
/// Hence, does not have 2D and 3D child classes.
/// Later, will be defined as a reference object type.
/// See ReferenceObject header file docstring for details.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Line_HeaderFile
#define _Line_HeaderFile

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
	class Line : public GeometryObject
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
		friend class Circle;
		friend class Plane;

	public:
		// ctor / dtor / operators
		Line(
			const std::shared_ptr<Point3D>& theEndPoint0,
			const std::shared_ptr<Point3D>& theEndPoint1);

		Line(const Line& rhs) = default;
		Line& operator=(const Line& rhs) = default;
		Line(Line&& rhs) noexcept = default;
		Line& operator=(Line&& rhs) noexcept = default;
		~Line() final = default;

		bool operator==(const Line& rhs) const;
		bool operator!=(const Line& rhs) const;
		bool operator+=(const Line& rhs) const;
		bool operator-=(const Line& rhs) const;

		// Methods:
		bool is2D() const;
		bool is3D() const;
		bool equalsGeometrically(ARGCOPY(Line) theLine) const;

		bool includes(ARGCOPY(Point3D) thePoint) const;
		bool intersects(ARGCOPY(Axis) theAxis) const;
		bool intersects(ARGCOPY(Line) theLine) const;
		bool coincides(ARGCOPY(Axis) theAxis) const;
		bool coincides(ARGCOPY(Line) theLine) const;
		bool isSkew(ARGCOPY(Axis) theAxis) const;
		bool isSkew(ARGCOPY(Line) theLine) const;
		auto intersect(ARGCOPY(Axis) theAxis) const;
		auto intersect(ARGCOPY(Line) theLine) const -> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>;
		auto project(ARGCOPY(Point3D) thePoint) const -> std::shared_ptr<Point3D>;
		double calculateDistance(ARGCOPY(Point3D) thePoint) const;
		double calculateDistance(ARGCOPY(Axis) theAxis) const;
		double calculateDistance(ARGCOPY(Line) theLine) const;

		auto getEndPoint0() const -> std::shared_ptr<Point3D>;
		auto getEndPoint1() const -> std::shared_ptr<Point3D>;
		auto getAxis() const->std::shared_ptr<Axis>;
		auto getDirectionVector() const->std::shared_ptr<Vector3D>;
		double getLength() const;
		auto getCommonReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>;
		void setEndPoint0(const std::shared_ptr<Point3D>& theEndPoint0);
		void setEndPoint1(const std::shared_ptr<Point3D>& theEndPoint1);
		auto createMidpoint() const -> std::shared_ptr<Point3D>;

	private:
		// Private default ctor used for cloning the object
		Line() = default;

		// Private methods
		void inspectEndPoints(
			const std::shared_ptr<Point3D>& theEndPoint0,
			const std::shared_ptr<Point3D>& theEndPoint1);

		// Members
		std::shared_ptr<Point3D> c_endPoint0 = nullptr;
		std::shared_ptr<Point3D> c_endPoint1 = nullptr;
	};
}

#endif
