/// <summary>
/// CoordSystem defines a coordinate system (CS) with an origin point and vectors for each of x, y and z axes.
/// 
/// INVARIANT:
///   No rule for the origim point coordinates.
///   Right-hand-rule is applicable to the axes vectors.
///   This is achieved by the constructors
///   even if the inputs (points or vectors) do not satisfy the rule.
/// 
/// For the origin point and axes,
/// point coords and vector components are prefered
/// rather than having a Point and Vector members
/// in order to prevent circular dependency.
/// Because Point and Vector objects (inheritting ReferenceObject class)
/// are defined with a CoordSystem member as the reference CS.
/// Hence, having Point and Vectors as members would create a circular dependency.
/// 
/// Point coords and axis vector components are wrt the global CS.
/// 
/// A 3D object currently.
/// Hence, does not have 2D and 3D child classes.
/// 
/// See docstring of GlobalCoordSystem.hxx for the global CS.
/// 
/// See GeometryObject.hxx for the details about this library.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _CoordSystem_HeaderFile
#define _CoordSystem_HeaderFile

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
	class CoordSystem : public GeometryObject
	{
		friend class GeometryObject;
		friend class ReferenceObject;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

		// Members
		std::array<double, 3> c_originCoords = { {} };
		std::array<double, 3> c_axisComponentsX = { {} };
		std::array<double, 3> c_axisComponentsY = { {} };
		std::array<double, 3> c_axisComponentsZ = { {} };
		bool c_isGlobal = false;

	public:
		// ctor / dtor / operators
		CoordSystem(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Point3D) thePointOnAxisX,
			ARGCOPY(Point3D) thePointOnAxisY);
		CoordSystem(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Vector3D) theAxisVectorX,
			ARGCOPY(Vector3D) theAxisVectorY);
		CoordSystem(
			const std::array<double, 3>& theOriginCoords,
			const std::array<double, 3>& theAxisComponentsX,
			const std::array<double, 3>& theAxisComponentsY,
			const std::array<double, 3>& theAxisComponentsZ);

		CoordSystem(const CoordSystem& rhs) = default;
		CoordSystem& operator=(const CoordSystem& rhs) = default;
		CoordSystem(CoordSystem&& rhs) noexcept = default;
		CoordSystem& operator=(CoordSystem&& rhs) noexcept = default;
		~CoordSystem() override = default;

		bool operator==(const CoordSystem& rhs) const;
		bool operator!=(const CoordSystem& rhs) const;
		bool operator+=(const CoordSystem& rhs) const;
		bool operator-=(const CoordSystem& rhs) const;

		// Methods
		bool is2D() const;
		bool is3D() const;
		bool equalsGeometrically(ARGCOPY(CoordSystem) theCoordSystem) const;

		auto getOriginCoords() const -> std::array<double, 3>;
		auto getAxisComponentsX() const -> std::array<double, 3>;
		auto getAxisComponentsY() const -> std::array<double, 3>;
		auto getAxisComponentsZ() const -> std::array<double, 3>;
		auto getAxisAsVectorX() const -> std::shared_ptr<Vector3D>;
		auto getAxisAsVectorY() const -> std::shared_ptr<Vector3D>;
		auto getAxisAsVectorZ() const -> std::shared_ptr<Vector3D>;
		auto getAxesAsVector() const->std::vector<std::shared_ptr<Vector3D>>;
		inline auto getIsGlobal() const { return c_isGlobal; };
		void setOriginCoords(const std::array<double, 3>& theOriginCoords);

		bool isGlobal() const;
		bool isParallel(ARGCOPY(CoordSystem) theCoordSystem) const;
		bool isIdentical(ARGCOPY(CoordSystem) theCoordSystem) const;

		auto createPoint(const std::array<double, 3>& theCoords);
		auto measurePointCoords(ARGCOPY(PointBase) thePoint) const -> std::array<double, 3>;
		auto measureVectorComponents(ARGCOPY(VectorBase) theVector) const -> std::array<double, 3>;
		auto rotatePointAboutAxisX(ARGCOPY(PointBase) thePoint, double theAngle);
		auto rotatePointAboutAxisY(ARGCOPY(PointBase) thePoint, double theAngle);
		auto rotatePointAboutAxisZ(ARGCOPY(PointBase) thePoint, double theAngle);

	protected:
		inline void setIsGlobal(bool theIsGlobal) { c_isGlobal = theIsGlobal; };

	private:
		// Private default ctor used for cloning the object
		CoordSystem() = default;

		// Private methods
		void setMembers(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Point3D) thePointX,
			ARGCOPY(Point3D) thePointY);
		void setMembers(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Vector3D) theAxisVectorX,
			ARGCOPY(Vector3D) theAxisVectorY);
		void setMembers(
			const std::array<double, 3>& theOriginCoords,
			const std::array<double, 3>& theAxisVectorX,
			const std::array<double, 3>& theAxisVectorY,
			const std::array<double, 3>& theAxisZ);
		auto rotatePointBase(
			ARGCOPY(PointBase) thePoint,
			const double& theAngle,
			int theAxis0,
			int theAxis1,
			int theAxis2) -> std::shared_ptr<Point3D>;
	};
}

#endif
