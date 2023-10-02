/// <summary>
/// The implementation base class for the the classes defined wrt a reference coordinate system (CS).
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
/// Point and Vector are the two reference object types.
/// 
/// Reference objects have 2D and 3D types (e.g. Point2D and Vector3D).
/// Other library objects do not have 2D and 3D types
/// as a reference CS is required to define the dimensionality.
/// 
/// The z dimension (i.e. coord for a point and component for a vector) is zero for 2D objects by default.
/// 
/// Axis, Line, Circle and Plane classes inherit from the base class (GeometryObject)
/// and do not have 2D and 3D types (e.g. Axis2D) currently.
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
/// Hence, currently, the types with more than one ReferenceObject members (Axis, Line, Circle and Plane)
/// do not have a reference CS and 2D/3D types
/// and inherit from the root class (GeometryObject).
/// Instead, they have three member functions:
/// is2D, is3D and getCommonReferenceCoordSystem.
/// 
/// LATER, setReferenceCoordSystem METHOD OF ReferenceObject WILL BE REMOVED
/// AND AXIS, LINE AND CIRCLE CLASSES WILL BE CHILD CLASS OF ReferenceObject.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#pragma warning(disable : 4290)

#ifndef _ReferenceObject_HeaderFile
#define _ReferenceObject_HeaderFile

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
	class ReferenceObject : public GeometryObject
	{
		friend class GeometryObject;
		friend class CoordSystem;
		friend class GlobalCoordSystem;
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
		int c_dimensionCount = GeometryParameters::DIMENSIONS::D3;
		std::shared_ptr<CoordSystem> c_referenceCoordSystem = nullptr;

		// ctor / dtor / operators
		explicit ReferenceObject(const int theDimensionCount);
		explicit ReferenceObject(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
		ReferenceObject(
			const int theDimensionCount,
			const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);

	public:
		ReferenceObject() = default;
		ReferenceObject(const ReferenceObject& rhs) = default;
		ReferenceObject& operator=(const ReferenceObject& rhs) = default;
		ReferenceObject(ReferenceObject&& rhs) = default;
		ReferenceObject& operator=(ReferenceObject&& rhs) = default;
		~ReferenceObject() override = default;

		bool operator==(const ReferenceObject& rhs) const;
		bool operator!=(const ReferenceObject& rhs) const;
		bool operator+=(const ReferenceObject& rhs) const;
		bool operator-=(const ReferenceObject& rhs) const;

		// Methods
		bool is2D() const;
		bool is3D() const;
		static bool is2DStrict(
			ARGCOPY(ReferenceObject) theReferenceObject0,
			ARGCOPY(ReferenceObject) theReferenceObject1);

		bool equalsRef(ARGCOPY(ReferenceObject) theReferenceObject) const;
		bool equalsGeometricallyRef(ARGCOPY(ReferenceObject) theReferenceObject) const;

		auto getReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>;
		auto getPointWithMyCoordSystem(ARGCOPY(PointBase) thePoint) const -> std::shared_ptr<Point3D>;
		auto getVectorWithMyCoordSystem(ARGCOPY(VectorBase) theVector) const -> std::shared_ptr<Vector3D>;

	protected:
		int getDimensionCount() const;
		void setDimensionCount(const int theDimensionCount);
		void setMembersRef(const int theDimensionCount, const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
		void setReferenceCoordSystemBase(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem);
	};
}

#endif
