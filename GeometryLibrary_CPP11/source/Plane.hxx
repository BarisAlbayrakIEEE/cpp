/// <summary>
/// Plane defines a plane with a passing point and normal vector.
/// See GeometryObject.hxx for the details about this library
/// 
/// INVARIANT:
///   Cluster of infinite number of points created by translating a point by all vectors pepandicular to a vector.
///   The point is called the passing point,
///   and the vector is called the normal vector.
/// 
/// A 3D object currently.
/// Hence, does not have 2D and 3D child classes.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _Plane_HeaderFile
#define _Plane_HeaderFile

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
	class Plane : public GeometryObject
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
		friend class Axis;
		friend class Line;
		friend class Circle;

	public:
		// ctor / dtor / operators
		explicit Plane(const std::array<double, 4>& theEC);
		explicit Plane(const std::vector<double, std::allocator<double>>& theEC);
		Plane(
			const std::shared_ptr<Point3D>& thePassingPoint,
			const std::shared_ptr<Vector3D>& theNormalVector);
		Plane(
			const std::shared_ptr<Point3D>& thePassingPoint,
			ARGCOPY(Vector3D) theInPlaneVectorO,
			ARGCOPY(Vector3D) theInPlaneVector1);
		Plane(
			const std::shared_ptr<Point3D>& thePassingPoint,
			ARGCOPY(Point3D) thePoint1,
			ARGCOPY(Point3D) thePoint2);

		Plane(const Plane& rhs) = default;
		Plane& operator=(const Plane& rhs) = default;
		Plane(Plane&& rhs) noexcept = default;
		Plane& operator=(Plane&& rhs) noexcept = default;
		~Plane() final = default;

		bool operator==(const Plane& rhs) const;
		bool operator!=(const Plane& rhs) const;
		bool operator+=(const Plane& rhs) const;
		bool operator-=(const Plane& rhs) const;

		// Methods
		bool is2D() const;
		bool is3D() const;
		bool equalsGeometrically(ARGCOPY(Plane) thePlane) const;

		auto getPassingPoint() const -> std::shared_ptr<Point3D>;
		auto getNormalVector() const -> std::shared_ptr<Vector3D>;
		auto getEC() const -> std::array<double, 4>;
		auto getCommonReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>;
		void setMembers(const std::array<double, 4>& theEC);
		void setMembers(const std::vector<double, std::allocator<double>>& theEC);
		void setMembers(
			const std::shared_ptr<Point3D>& thePassingPoint,
			const std::shared_ptr<Vector3D>& theNormalVector);
		void setPassingPoint(const std::shared_ptr<Point3D>& thePassingPoint);
		void setNormalVector(const std::shared_ptr<Vector3D>& theNormalVector);

		bool intersects(ARGCOPY(Axis) theAxis) const;
		bool intersects(ARGCOPY(Line) theLine) const;
		bool includes(ARGCOPY(Point3D) thePoint) const;
		bool includes(ARGCOPY(Axis) theAxis) const;
		bool includes(ARGCOPY(Line) theLine) const;
		auto calculateCoordX(const double& theCoordY, const double& theCoordZ) const -> std::pair<bool, double>;
		auto calculateCoordY(const double& theCoordZ, const double& theCoordX) const -> std::pair<bool, double>;
		auto calculateCoordZ(const double& theCoordX, const double& theCoordY) const -> std::pair<bool, double>;
		auto calculateCoordXY(const double& theCoordZ) const;
		auto calculateCoordYZ(const double& theCoordX) const;
		auto calculateCoordZX(const double& theCoordY) const;
		double calculateDistance(ARGCOPY(Plane) thePlane) const;
		double calculateDistance(ARGCOPY(Point3D) thePoint) const;
		double calculateDistance(ARGCOPY(Axis) theAxis) const;
		double calculateDistance(ARGCOPY(Line) theLine) const;
		auto intersect(ARGCOPY(Plane) thePlane) const -> std::pair<GeometryParameters::INTERSECTION2, std::shared_ptr<Axis>>;
		auto intersect(ARGCOPY(Axis) theAxis) const -> std::pair<GeometryParameters::INTERSECTION2, std::shared_ptr<Point3D>>;
		auto intersect(ARGCOPY(Line) theLine) const -> std::pair<GeometryParameters::INTERSECTION2, std::shared_ptr<Point3D>>;
		auto project(ARGCOPY(Point3D) thePoint) const -> std::shared_ptr<Point3D>;
		auto project(ARGCOPY(Vector3D) theVector) const -> std::shared_ptr<Vector3D>;
		auto project(ARGCOPY(Axis) theAxis) const -> std::shared_ptr<Axis>;
		auto project(ARGCOPY(Line) theLine) const -> std::shared_ptr<Line>;
		auto createPoint(const double& theFactor) const;

	private:
		// Private default ctor used for cloning the object
		Plane() = default;

		// Private methods
		void inspectEC(const std::array<double, 4>& theEC) const;
		void inspectEC(const std::vector<double, std::allocator<double>>& theEC) const;
		void applyEC(const std::array<double, 4>& theEC);
		void applyEC(const std::vector<double, std::allocator<double>>& theEC);

		// Members
		std::shared_ptr<Point3D> c_passingPoint = nullptr;
		std::shared_ptr<Vector3D> c_normalVector = nullptr;
	};
}

#endif
