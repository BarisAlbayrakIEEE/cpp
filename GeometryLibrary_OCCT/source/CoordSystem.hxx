/// <summary>
/// CoordSystem defines a coordinate system (CS) with an origin point and vectors for each of x, y and z axes.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
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

#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Handle_HeaderFile
#include <Standard_Handle.hxx>
#endif
#ifndef _Standard_Type_HeaderFile
#include <Standard_Type.hxx>
#endif
#ifndef _Standard_Size_HeaderFile
#include <Standard_Size.hxx>
#endif
#ifndef _GeometryAbstractObject_HeaderFile
#include "GeometryAbstractObject.hxx"
#endif
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

	DEFINE_STANDARD_HANDLE(CoordSystem, GeometryObject)

	class CoordSystem : public virtual GeometryAbstractObject, public GeometryObject
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
		/*
		 * DEF: Modifications on OCCT
		*/
		void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(CoordSystem);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT bool is2D() const;
		Standard_EXPORT bool is3D() const;
		Standard_EXPORT virtual void Destroy();

	public:
		Standard_EXPORT bool equals(ARGCOPY(CoordSystem) theCoordSystem) const;
		Standard_EXPORT bool equalsGeometrically(ARGCOPY(CoordSystem) theCoordSystem) const;

		// Members
	protected:
		bool c_isGlobal;

	private:
		arrayS3 c_originCoords;
		arrayS3 c_axisComponentsX;
		arrayS3 c_axisComponentsY;
		arrayS3 c_axisComponentsZ;

	public:
		Standard_EXPORT CoordSystem(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Point3D) thePointOnAxisX,
			ARGCOPY(Point3D) thePointOnAxisY);
		Standard_EXPORT CoordSystem(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Vector3D) theAxisVectorX,
			ARGCOPY(Vector3D) theAxisVectorY);

	protected:
		Standard_EXPORT CoordSystem(
			const arrayS3& theOriginCoords,
			const arrayS3& theAxisComponentsX,
			const arrayS3& theAxisComponentsY,
			const arrayS3& theAxisComponentsZ);

	public:
		// ctor / dtor / operators
		Standard_EXPORT CoordSystem(const CoordSystem& rhs);
		Standard_EXPORT CoordSystem& operator=(const CoordSystem& rhs);
		Standard_EXPORT CoordSystem(CoordSystem&& rhs) noexcept;
		Standard_EXPORT CoordSystem& operator=(CoordSystem&& rhs) noexcept;
		Standard_EXPORT bool operator==(const CoordSystem& rhs);
		Standard_EXPORT bool operator!=(const CoordSystem& rhs);
		Standard_EXPORT bool operator+=(const CoordSystem& rhs);
		Standard_EXPORT bool operator-=(const CoordSystem& rhs);
		Standard_EXPORT ~CoordSystem();

	private:
		Standard_EXPORT void copyBase(const CoordSystem& rhs);

	public:
		DEFINE_STANDARD_RTTIEXT(CoordSystem, GeometryObject)

		// Methods
	public:
		Standard_EXPORT arrayS3 getOriginCoords() const;
		Standard_EXPORT arrayS3 getAxisComponentsX() const;
		Standard_EXPORT arrayS3 getAxisComponentsY() const;
		Standard_EXPORT arrayS3 getAxisComponentsZ() const;
		Standard_EXPORT OUTVAL(Vector3D) getAxisAsVectorX() const;
		Standard_EXPORT OUTVAL(Vector3D) getAxisAsVectorY() const;
		Standard_EXPORT OUTVAL(Vector3D) getAxisAsVectorZ() const;
		Standard_EXPORT std::vector<Handle(Vector3D)> getAxesAsVector() const;
		Standard_EXPORT void setOriginCoords(const arrayS3& theOriginCoords);

		Standard_EXPORT bool isGlobal() const;
		Standard_EXPORT bool isParallel(ARGCOPY(CoordSystem) theCoordSystem) const;
		Standard_EXPORT bool isIdentical(ARGCOPY(CoordSystem) theCoordSystem) const;

		Standard_EXPORT OUTVAL(PointBase) createPoint(int theDimenslonCount, const arrayS3& theCoords);
		Standard_EXPORT arrayS3 measurePointCoords(ARGCOPY(PointBase) thePoint) const;
		Standard_EXPORT arrayS3 measureVectorComponents(ARGCOPY(VectorBase) theVector) const;
		Standard_EXPORT OUTVAL(Point3D) rotatePointAboutAxisX(ARGCOPY(PointBase) thePoint, double theAngle);
		Standard_EXPORT OUTVAL(Point3D) rotatePointAboutAxisY(ARGCOPY(PointBase) thePoint, double theAngle);
		Standard_EXPORT OUTVAL(Point3D) rotatePointAboutAxisZ(ARGCOPY(PointBase) thePoint, double theAngle);

	private:
		Standard_EXPORT void setMembers(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Point3D) thePointX,
			ARGCOPY(Point3D) thePointY);
		Standard_EXPORT void setMembers(
			ARGCOPY(Point3D) theOriginPoint,
			ARGCOPY(Vector3D) theAxisVectorX,
			ARGCOPY(Vector3D) theAxisVectorY);
		Standard_EXPORT void setMembers(
			const arrayS3& theOriginCoords,
			const arrayS3& theAxisVectorX,
			const arrayS3& theAxisVectorY,
			const arrayS3& theAxisZ);
		Standard_EXPORT OUTVAL(Point3D) rotatePointBase(
			ARGCOPY(PointBase) thePoint,
			const double& theAngle,
			int theAxisO,
			int theAxis1,
			int theAxis2);
	};
}

#endif
