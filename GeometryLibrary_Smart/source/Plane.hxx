/// <summary>
/// Plane defines a plane with a passing point and normal vector.
/// See GeometryObject.hxx for details
/// 
/// NO INVARIANT
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
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

	class Plane : public GeometryObject
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
		friend class Line;
		friend class Circle;

	public:
		bool is2D() const;
		bool is3D() const;

	public:
		bool equals(ARGCOPY(Plane) thePlane) const;
		bool equalsGeometrically(ARGCOPY(Plane) thePlane) const;

	private:
		// Members
		std::shared_ptr<PointBase> c_passingPoint;
		std::shared_ptr<VectorBase> c_normalVector;
		arrayS4 c_EC; // C0x + C1y + C2z + C3 = O

	public:
		// ctor / dtor / operators
		Plane(const arrayS4& theEC);
		Plane(const vectorInput1D& theEC);
		Plane(
			PointBase& thePoint,
			VectorBase& theNormalVector);
		Plane(
			PointBase& thePoint,
			ARGCOPY(VectorBase) theInPlaneVectorO,
			ARGCOPY(VectorBase) theInPlaneVector1);
		Plane(
			PointBase& thePoint0,
			ARGCOPY(PointBase) thePoint1,
			ARGCOPY(PointBase) thePoint2);

		Plane(const Plane& rhs);
		Plane& operator=(const Plane& rhs);
		Plane(Plane&& rhs) noexcept;
		Plane& operator=(Plane&& rhs) noexcept;
		bool operator==(const Plane& rhs);
		bool operator!=(const Plane& rhs);
		bool operator+=(const Plane& rhs);
		bool operator-=(const Plane& rhs);
		~Plane();

	private:
		void copyBase(const Plane& rhs);
		void Destroy();

			// Methods
	public:
		PointBase& getPassingPoint() const;
		VectorBase& getNormalVector() const;
		arrayS4 getEC() const;
		CoordSystem* getCommonReferenceCoordSystem() const;
		void setMembers(const arrayS4& theEC);
		void setMembers(const vectorInput1D& theEC);
		void setMembers(
			PointBase& thePassingPoint,
			VectorBase& theNormalVector);
		void setPassingPoint(PointBase& thePassingPoint);
		void setNormalVector(VectorBase& theNormalVector);
		void setEC(const arrayS4& theEC);
		void setEC(const vectorInput1D& theEC);

		bool intersects(ARGCOPY(Axis) theAxis) const;
		bool intersects(ARGCOPY(Line) theLine) const;
		bool includes(ARGCOPY(PointBase) thePoint) const;
		bool includes(ARGCOPY(Axis) theAxis) const;
		bool includes(ARGCOPY(Line) theLine) const;
		std::pair<bool, double> calculateCoordX(const double& theCoordY, const double& theCoordZ) const;
		std::pair<bool, double> calculateCoordY(const double& theCoordZ, const double& theCoordX) const;
		std::pair<bool, double> calculateCoordZ(const double& theCoordX, const double& theCoordY) const;
		std::pair<bool, std::array<double, 2>> calculateCoordXY(const double& theCoordZ) const;
		std::pair<bool, std::array<double, 2>> calculateCoordYZ(const double& theCoordX) const;
		std::pair<bool, std::array<double, 2>> calculateCoordZX(const double& theCoordY) const;
		double calculateDistance(ARGCOPY(Plane) thePlane) const;
		double calculateDistance(ARGCOPY(PointBase) thePoint) const;
		double calculateDistance(ARGCOPY(Axis) theAxis) const;
		double calculateDistance(ARGCOPY(Line) theLine) const;
		std::pair<INTERSECTION2, Axis*> intersect(ARGCOPY(Plane) thePlane) const;
		std::pair<INTERSECTION2, Point3D*> intersect(ARGCOPY(Axis) theAxis) const;
		std::pair<INTERSECTION2, Point3D*> intersect(ARGCOPY(Line) theLine) const;
		Point3D project(ARGCOPY(PointBase) thePoint) const;
		Vector3D* project(ARGCOPY(VectorBase) theVector) const;
		Axis* project(ARGCOPY(Axis) theAxis) const;
		Line* project(ARGCOPY(Line) theLine) const;
		Point3D createPoint(const double& theFactor) const;

	private:
		bool equalsBase(ARGCOPY(Plane) thePlane) const;
		void inspectEC(const arrayS4& theEC);
		void inspectEC(const vectorInput1D& theEC);
		void updateEC();
		void applyEC();
	};
}

#endif
