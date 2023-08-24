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
/// A 3D object by default.
/// Hence, does not have 2D and 3D child classes.
/// See GeometryObject.hxx for the details.
/// 
/// See docstring of GlobalCoordSystem.hxx for the global CS.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
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

	class CoordSystem : public GeometryObject
	{
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

	public:
		bool is2D() const;
		bool is3D() const;

	public:
		bool equals(ARGCOPY(CoordSystem) theCoordSystem) const;
		bool equalsGeometrically(ARGCOPY(CoordSystem) theCoordSystem) const;

		// Members
	protected:
		bool c_isGlobal;

	private:
		arrayS3 c_originCoords;
		arrayS3 c_axisComponentsX;
		arrayS3 c_axisComponentsY;
		arrayS3 c_axisComponentsZ;

	public:
		CoordSystem(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Point3D) thePointOnAxisX,
			ARGCOPY(Point3D) thePointOnAxisY);
		CoordSystem(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Vector3D) theAxisVectorX,
			ARGCOPY(Vector3D) theAxisVectorY);

	protected:
		CoordSystem(
			const arrayS3& theOriginCoords,
			const arrayS3& theAxisComponentsX,
			const arrayS3& theAxisComponentsY,
			const arrayS3& theAxisComponentsZ);

	public:
		// ctor / dtor / operators
		CoordSystem(const CoordSystem& rhs);
		CoordSystem& operator=(const CoordSystem& rhs);
		CoordSystem(CoordSystem&& rhs) noexcept;
		CoordSystem& operator=(CoordSystem&& rhs) noexcept;
		bool operator==(const CoordSystem& rhs);
		bool operator!=(const CoordSystem& rhs);
		bool operator+=(const CoordSystem& rhs);
		bool operator-=(const CoordSystem& rhs);
		~CoordSystem();

	private:
		void copyBase(const CoordSystem& rhs);
		void Destroy();

		// Methods
	public:
		arrayS3 getOriginCoords() const;
		arrayS3 getAxisComponentsX() const;
		arrayS3 getAxisComponentsY() const;
		arrayS3 getAxisComponentsZ() const;
		Vector3D getAxisAsVectorX() const;
		Vector3D getAxisAsVectorY() const;
		Vector3D getAxisAsVectorZ() const;
		std::vector<Vector3D> getAxesAsVector() const;
		void setOriginCoords(const arrayS3& theOriginCoords);

		bool isGlobal() const;
		bool isParallel(ARGCOPY(CoordSystem) theCoordSystem) const;
		bool isIdentical(ARGCOPY(CoordSystem) theCoordSystem) const;

		PointBase createPoint(int theDimenslonCount, const arrayS3& theCoords);
		arrayS3 measurePointCoords(ARGCOPY(PointBase) thePoint) const;
		arrayS3 measureVectorComponents(ARGCOPY(VectorBase) theVector) const;
		Point3D rotatePointAboutAxisX(ARGCOPY(PointBase) thePoint, double theAngle);
		Point3D rotatePointAboutAxisY(ARGCOPY(PointBase) thePoint, double theAngle);
		Point3D rotatePointAboutAxisZ(ARGCOPY(PointBase) thePoint, double theAngle);

	private:
		void setMembers(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Point3D) thePointX,
			ARGCOPY(Point3D) thePointY);
		void setMembers(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Vector3D) theAxisVectorX,
			ARGCOPY(Vector3D) theAxisVectorY);
		void setMembers(
			const arrayS3& theOriginCoords,
			const arrayS3& theAxisVectorX,
			const arrayS3& theAxisVectorY,
			const arrayS3& theAxisZ);
		Point3D rotatePointBase(
			ARGCOPY(PointBase) thePoint,
			const double& theAngle,
			int theAxisO,
			int theAxis1,
			int theAxis2);
	};
}

#endif
