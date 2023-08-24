/// <summary>
/// Line defines a line piece between two points.
/// Cluster of the points between two end points.
/// 
/// INVARIANT:
///   The end points shall be geometrically unequal (i.e. non-coincident).
/// 
/// The points can be defined wrt different reference CSs.
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

	class Line : public GeometryObject
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
		bool is2D() const;
		bool is3D() const;
		void Destroy();

	public:
		bool equals(ARGCOPY(Line) theLine) const;
		bool equalsGeometrically(ARGCOPY(Line) theLine) const;

	private:
		// Members
		std::shared_ptr<Axis> c_axis;
		std::shared_ptr<PointBase> c_endPoint0;
		std::shared_ptr<PointBase> c_endPoint1;
		double c_length;

	public:
		// ctor / dtor / operators
		Line(
			PointBase& theEndPoint0,
			PointBase& theEndPoint1);

		Line(const Line& rhs);
		Line& operator=(const Line& rhs);
		Line(Line&& rhs) noexcept;
		Line& operator=(Line&& rhs) noexcept;
		bool operator==(const Line& rhs);
		bool operator!=(const Line& rhs);
		bool operator+=(const Line& rhs);
		bool operator-=(const Line& rhs);
		~Line();

	private:
		void copyBase(const Line& rhs);

		// Methods:
	public:
		bool includes(ARGCOPY(PointBase) thePoint) const;
		bool intersects(ARGCOPY(Axis) theAxis) const;
		bool intersects(ARGCOPY(Line) theLine) const;
		bool coincides(ARGCOPY(Axis) theAxis) const;
		bool coincides(ARGCOPY(Line) theLine) const;
		bool isSkew(ARGCOPY(Axis) theAxis) const;
		bool isSkew(ARGCOPY(Line) theLine) const;
		std::pair<INTERSECTION1, PointBase*> intersect(ARGCOPY(Axis) theAxis) const;
		std::pair<INTERSECTION1, PointBase*> intersect(ARGCOPY(Line) theLine) const;
		PointBase project(ARGCOPY(PointBase) thePoint);
		double calculateDistance(ARGCOPY(PointBase) thePoint);
		double calculateDistance(ARGCOPY(Axis) theAxis);
		double calculateDistance(ARGCOPY(Line) theLine);

		Axis& getAxis() const;
		VectorBase& getDirectionVector() const;
		std::vector<PointBase> getEndPoints() const;
		PointBase& getEndPoint0() const;
		PointBase& getEndPoint1() const;
		double getLength() const;
		CoordSystem* getCommonReferenceCoordSystem() const;
		void setAxis(Axis& theAxis);
		void setEndPoint0(PointBase& theEndPoint0);
		void setEndPoint1(PointBase& theEndPoint1);
		PointBase createMidpoint() const;

	private:
		bool equalsBase(ARGCOPY(Line) theLine) const;
	};
}

#endif
