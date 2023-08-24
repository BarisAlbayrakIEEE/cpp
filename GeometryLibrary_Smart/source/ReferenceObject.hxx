/// <summary>
/// The implementation base class for the the classes defined wrttt a reference coordinate system (CS).
/// 
/// This project is a simple implementation for the geometry objects (e.g. point, vector and circle)
/// and interactions of these objects.
/// Complex geometries (e.g. a geometrically undefined surface) are not included yet.
/// 
/// The reference objects are defined wrt to a CS in order to
///   1. simplify the tracebility of the objects
///   2. switch between spaces, especially 3D to 2D (or 2D to 3D)
/// The 1st functionality is a fundamental issue in engineering applications.
/// For example, an assembly of bodies is defined in a global CS
/// while each body has its own CS.
/// The 2nd functionality is similar to the 1st one.
/// Its benefit is more clear when it comes to switching between 2D and 3D spaces.
/// All CAE software (e.g. Catia) supports switching to 2D space.
/// This library can be utilized in this respect.
/// 
/// Currently, Point and Vector are the two reference object types.
/// 
/// Reference objects have 2D and 3D types (e.g. Point2D and Vector3D).
/// Other library objects do not have 2D and 3D types
/// as a reference CS is required to define the dimensionality.
/// 
/// The z dimension (i.e. coord for a point and component for a vector) is zero for 2D objects by defaullt.
/// 
/// Axis, Line, Circle and Plane classes inherit from the base class (GeometryObject)
/// and do not have 2D and 3D types (e.g. Axis2D).
/// Actually, only the classes inheritting from ReferenceObject have 2D and 3D types.
/// The problem is not 2D/3D types but the reference CS.
/// An Axis cannot have a reference CS
/// because it has two members (point and vector)
/// and the reference CSs of the two may not be the same.
/// This can be satisfied as a precondition in the ctors of Axis
/// but the reference CSs of the point and vector members can be modified
/// after the Axis object instantiation.
/// Hence, setReferenceCoordSystem method of ReferenceObject shall be removed
/// in order for Axis to have a reference CS
/// But this solution is not preferable.
/// Hence, the types with more than one ReferenceObject members (Axis, Line, Circle and Plane)
/// do not have a reference CS and 2D/3D types
/// and inherit from the root class (GeometryObject).
/// Instead, they have three member functions:
/// is2D, is3D and getCommonReferenceCoordSystem
/// 
/// The objects in the library are defined in an implementation hyerarchy.
/// The slicing problem that comes with the inheritance is
/// solved by deleting the copy/move ctors and operators.
/// Hence, copy/move ctors and aassignments of all base classes accept for PointBase
/// (GeometryObject, ReferenceObject, CoordSystem, VectorBase) are deleted.
/// Point2D and Point3D classes do not have additional members
/// so that the slicing issue is not a problem for them.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _ReferenceObject_HeaderFile
#define _ReferenceObject_HeaderFile

#ifndef _ReferenceAbstractObject_HeaderFile
#include "ReferenceAbstractObject.hxx"
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

	class ReferenceObject : public virtual ReferenceAbstractObject, public GeometryObject
	{
	public:
		bool is2D() const;
		bool is3D() const;
		static bool is2DStrict(
			ARGCOPY(ReferenceObject) theReferenceObject0,
			ARGCOPY(ReferenceObject) theReferenceObject1);

	public:
		virtual bool equals(ARGCOPY(ReferenceObject) theReferenceObject) const;
		virtual bool equalsGeometrically(ARGCOPY(ReferenceObject) theReferenceObject) const;

	protected:
		// Members
		int c_dimensionCount;
		std::shared_ptr<CoordSystem> c_referenceCoordSystem;

		// ctor / dtor / operators
	protected:
		ReferenceObject(const int theDimensionCount);
		ReferenceObject(CoordSystem& theReferenceCoordSystem);
		ReferenceObject(
			const int theDimensionCount,
			CoordSystem& theReferenceCoordSystem);

	public:
		ReferenceObject(const ReferenceObject& rhs) = delete;
		ReferenceObject& operator=(const ReferenceObject& rhs) = delete;
		ReferenceObject(ReferenceObject&& rhs) = delete;
		ReferenceObject& operator=(ReferenceObject&& rhs) = delete;
		bool operator==(const ReferenceObject& rhs);
		bool operator!=(const ReferenceObject& rhs);
		bool operator+=(const ReferenceObject& rhs);
		bool operator-=(const ReferenceObject& rhs);
		virtual ~ReferenceObject();

	protected:
		virtual void copyBase(const ReferenceObject& rhs);
		void Destroy();

		// Methods
	public:
		CoordSystem& getReferenceCoordSystem() const;
		PointBase getPointWithMyCoordSystem(ARGCOPY(PointBase) thePoint) const;
		VectorBase getVectorWithMyCoordSystem(ARGCOPY(VectorBase) theVector) const;

	protected:
		void inspectDimensionCount(const int theDimensionCount) const;
		int getDimensionCount() const;
		void setDimensionCount(const int theDimensionCount);
		void setMembers(const int theDimensionCount, CoordSystem& theReferenceCoordSystem);
		void setReferenceCoordSystemBase(CoordSystem& theReferenceCoordSystem);
	};
}

#endif
