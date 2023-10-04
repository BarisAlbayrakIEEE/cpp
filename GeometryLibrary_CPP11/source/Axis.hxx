/// <summary>
/// Axis defines a directed infinite line.
/// 
/// INVARIANT:
///   Cluster of infinite number of points created by translating a point by a vector.
///   The translation vector is called the direction vector,
///   and the translated point is called the passing point.
/// 
/// Although an axis is infinite, it has a direction.
/// 
/// A 3D object currently.
/// Hence, does not have 2D and 3D child classes.
/// Later, will be defined as a reference object type.
/// See ReferenceObject header file docstring for details.
/// 
/// See GeometryObject.hxx for the details about this library.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Axis_HeaderFile
#define _Axis_HeaderFile

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
	class Axis : public GeometryObject
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class CoordSystem;
		friend class Point3D;
		friend class Point2D;
		friend class Point3D;
		friend class Vector3D;
		friend class Vector2D;
		friend class Vector3D;
		friend class Line;
		friend class Circle;
		friend class Plane;

	public:
		// ctor / dtor / operators
		Axis(
			std::shared_ptr<Point3D> thePassingPoint,
			std::shared_ptr<Vector3D> theDirectionVector);
		Axis(
			std::shared_ptr<Point3D> thePoint0,
			ARGCOPY(Point3D) thePoint1);
		explicit Axis(const std::array<std::array<double, 3>, 2>& theEC);
		explicit Axis(const std::vector<std::vector<double, std::allocator<double>>>& theEC);

		Axis(const Axis& rhs) = default;
		Axis& operator=(const Axis& rhs) = default;
		Axis(Axis&& rhs) noexcept = default;
		Axis& operator=(Axis&& rhs) noexcept = default;
		~Axis() final = default;

		bool operator==(const Axis& rhs) const;
		bool operator!=(const Axis& rhs) const;
		bool operator+=(const Axis& rhs) const;
		bool operator-=(const Axis& rhs) const;

		bool is2D() const;
		bool is3D() const;
		bool equalsGeometrically(ARGCOPY(Axis) theAxis) const;

		auto getPassingPoint() const -> std::shared_ptr<Point3D>;
		auto getDirectionVector() const -> std::shared_ptr<Vector3D>;
		auto getEC() const -> std::array<std::array<double, 3>, 2>;
		auto createPoint(const double& theFactor) const -> std::shared_ptr<Point3D>;
		auto getPassingPointAsVector() const;
		auto getPointAsVector(const double& theFactor) const;
		void setPassingPoint(const std::shared_ptr<Point3D>& thePassingPoint);
		void setDirectionVector(const std::shared_ptr<Vector3D>& theDirectionVector);
		bool isParallel(ARGCOPY(Axis) theAxis) const;
		bool isInTheSameDirection(ARGCOPY(Axis) theAxis) const;
		bool isNormal(ARGCOPY(Axis) theAxis) const;
		bool isSkew(ARGCOPY(Axis) theAxis) const;
		bool isSkew(ARGCOPY(Line) theLine) const;

		bool includes(ARGCOPY(Point3D) thePoint) const;
		bool intersects(ARGCOPY(Axis) theAxis) const;
		bool intersects(ARGCOPY(Line) theLine) const;
		bool coincides(ARGCOPY(Axis) theAxis) const;
		bool coincides(ARGCOPY(Line) theLine) const;
		auto intersect(ARGCOPY(Axis) theAxis) const-> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>;
		auto intersect(ARGCOPY(Line) theLine) const-> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>;
		auto project(ARGCOPY(Point3D) thePoint) const -> std::shared_ptr<Point3D>;
		double calculateDistance(ARGCOPY(Point3D) thePoint) const;
		double calculateDistance(ARGCOPY(Axis) theAxis) const;
		double calculateDistance(ARGCOPY(Line) theLine) const;
		auto findClosestPoints(ARGCOPY(Axis) theAxis) const -> std::vector<std::shared_ptr<Point3D>>;

		double calculateCoordX_fromCoordY(const double& theCoordY) const;
		double calculateCoordX_fromCoordZ(const double& theCoordZ) const;
		double calculateCoordY_fromCoordX(const double& theCoordX) const;
		double calculateCoordY_fromCoordZ(const double& theCoordZ) const;
		double calculateCoordZ_fromCoordX(const double& theCoordX) const;
		double calculateCoordZ_fromCoordY(const double& theCoordY) const;

	private:
		// Private default ctor used for cloning the object
		Axis() = default;

		// Private methods
		void setMembers(const std::array<std::array<double, 3>, 2>& theEC);
		void setMembers(const std::vector<std::vector<double, std::allocator<double>>>& theEC);
		void setMembers(
			const std::shared_ptr<Point3D>& thePassingPoint,
			const std::shared_ptr<Vector3D>& theDirectionVector);
		void inspectEC(const std::array<std::array<double, 3>, 2>& theEC) const;
		void inspectEC(const std::vector<std::vector<double, std::allocator<double>>>& theEC) const;
		void applyEC(const std::array<std::array<double, 3>, 2>& theEC);
		void applyEC(const std::vector<std::vector<double, std::allocator<double>>>& theEC);

		double calculateCoord(
			int theIndexCoordO,
			int theIndexCoordl,
			const double& theCoord) const;
		auto intersectBase(
			ARGCOPY(Axis) theAxis,
			const std::array<double, 3>& theCrossProductComponents) const
			-> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>;
		
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		// The axis Equation Coefficients(EC) of the form : (x - x0) / Vx = (y - y0) / Vy = (z - z0) / Vz
		// x0, y0 and z0 are coordinates of the passing point
		// Vx, Vy and Vz are the components of the direction vector
		// c_EC[0][i] = Passing point coord in ith direction
		// c_EC[1][i] = Direction vector component in ith direction
		std::shared_ptr<Point3D> c_passingPoint = nullptr;
		std::shared_ptr<Vector3D> c_directionVector = nullptr;
	};
}

#endif
