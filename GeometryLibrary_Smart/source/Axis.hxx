/// <summary>
/// Axis defines a directed infinite line.
/// Cluster of infinite number of points created by translating a point by a vector.
/// 
/// INVARIANT:
///   The translation vector is called the direction vector,
///   and the translated point is called the passing point.
/// 
/// Although an axis is infinite, it has a direction.
/// 
/// A 3D object by default.
/// Hence, does not have 2D and 3D child classes.
/// See GeometryObject.hxx for the details.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
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

	class Axis : public GeometryObject
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector2D;
		friend class Vector3D;
		friend class Line;
		friend class Circle;
		friend class Plane;

	public:
		bool is2D() const;
		bool is3D() const;

	public:
		bool equals(ARGCOPY(Axis) theAxis) const;
		bool equalsGeometrically(ARGCOPY(Axis) theAxis) const;

	private:
		// Members
		// The reference CS of the point member is assumed to be the reference CS.
		// The axis Equation Coefficients(EC) of the form : (x - x0) / Vx = (y - y0) / Vy = (z - z0) / Vz
		// x0, y0 and z0 are coordinates of the passing point
		// Vx, Vy and Vz are the components of the direction vector
		// c_EC[i][O] = Point coord in ith direction
		// c_EC[i][1] = Vector component in ith direction
		std::shared_ptr<PointBase> c_passingPoint;
		std::shared_ptr<VectorBase> c_directionVector;
		arrayS32 c_EC;

	public:
		// ctor / dtor / operators
		Axis(
			PointBase& thePassingPoint,
			VectorBase& theDirectionVector);
		Axis(
			PointBase& thePoint0,
			ARGCOPY(PointBase) thePoint1);
		Axis(
			Point2D& thePassingPoint,
			const double& theAngle);
		Axis(
			PointBase& thePassingPoint,
			const arrayS3& theAngles);
		Axis(const arrayS32& theEC);
		Axis(vectorInput2D theEC);

		Axis(const Axis& rhs);
		Axis& operator=(const Axis& rhs);
		Axis(Axis&& rhs) noexcept;
		Axis& operator=(Axis&& rhs) noexcept;
		bool operator==(const Axis& rhs);
		bool operator!=(const Axis& rhs);
		bool operator+=(const Axis& rhs);
		bool operator-=(const Axis& rhs);
		~Axis();

	private:
		void copyBase(const Axis& rhs);

	public:
		PointBase& getPassingPoint() const;
		VectorBase& getDirectionVector() const;
		arrayS32 getEC() const;
		PointBase createPoint(const double& theFactor) const;
		Vector3D getPassingPointAsVector() const;
		Vector3D getPointAsVector(const double& theFactor) const;
		CoordSystem* getCommonReferenceCoordSystem() const;
		void setPassingPoint(PointBase& thePassingPoint);
		void setDirectionVector(VectorBase& theDirectionVector);
		void setEC(const arrayS32& theEC);
		void setEC(const vectorInput2D& theEC);
		bool isParallel(ARGCOPY(Axis) theAxis) const;
		bool isParallel(
			ARGCOPY(Axis) theAxis,
			const double& theTolerance) const;
		bool isInTheSameDirection(ARGCOPY(Axis) theAxis) const;
		bool isInTheSameDirection(
			ARGCOPY(Axis) theAxis,
			const double& theTolerance) const;
		bool isNormal(ARGCOPY(Axis) theAxis) const;
		bool isNormal(
			ARGCOPY(Axis) theAxis,
			const double& theTolerance) const;
		bool isSkew(ARGCOPY(Axis) theAxis) const;
		bool isSkew(ARGCOPY(Line) theLine) const;

		bool includes(ARGCOPY(PointBase) thePoint) const;
		bool intersects(ARGCOPY(Axis) theAxis) const;
		bool intersects(ARGCOPY(Line) theLine) const;
		bool coincides(ARGCOPY(Axis) theAxis) const;
		bool coincides(ARGCOPY(Line) theLine) const;
		std::pair<INTERSECTION1, PointBase*> intersect(ARGCOPY(Axis) theAxis) const;
		std::pair<INTERSECTION1, PointBase*> intersect(ARGCOPY(Line) theLine) const;
		PointBase project(ARGCOPY(PointBase) thePoint) const;
		double calculateDistance(ARGCOPY(PointBase) thePoint) const;
		double calculateDistance(ARGCOPY(Axis) theAxis) const;
		double calculateDistance(ARGCOPY(Line) theLine) const;
		std::vector<PointBase> findClosestPoints(ARGCOPY(Axis) theAxis) const;

		double calculateCoordX_fromCoordY(const double& theCoordY) const;
		double calculateCoordX_fromCoordZ(const double& theCoordZ) const;
		double calculateCoordY_fromCoordX(const double& theCoordX) const;
		double calculateCoordY_fromCoordZ(const double& theCoordZ) const;
		double calculateCoordZ_fromCoordX(const double& theCoordX) const;
		double calculateCoordZ_fromCoordY(const double& theCoordY) const;

	private:
		void Destroy();
		void setMembers(const arrayS32& theEC);
		void setMembers(const vectorInput2D& theEC);
		void setMembers(
			PointBase& thePassingPoint,
			VectorBase& theDirectionVector);
		void inspectEC(const arrayS32& theEC);
		void inspectEC(const vectorInput2D& theEC);
		int determineDimensionCountUsingEC(const arrayS32& theEC) const;
		int determineDimensionCountUsingEC(const vectorInput2D& theEC) const;
		void updateEC();
		void applyEC();

		double calculateCoord(
			int theIndexCoordO,
			int theIndexCoordl,
			const double& theCoord) const;
		std::pair<INTERSECTION1, PointBase*> intersectBase(
			ARGCOPY(Axis) theAxis,
			const arrayS3& theCrossProductComponents) const;
		bool equalsBase(ARGCOPY(Axis) theAxis) const;
	};
}

#endif
