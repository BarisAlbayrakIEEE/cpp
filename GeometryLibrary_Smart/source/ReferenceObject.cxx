// baris.albayrak.ieee@gmail.com

#include "GeometryObject.hxx"
#include "ReferenceObject.hxx"
#include "CoordSystem.hxx"
#include "GlobalCoordSystem.hxx"
#include "PointBase.hxx"
#include "Point2D.hxx"
#include "Point3D.hxx"
#include "VectorBase.hxx"
#include "Vector2D.hxx"
#include "Vector3D.hxx"
#include "Axis.hxx"
#include "Line.hxx"
#include "Circle.hxx"
#include "Plane.hxx"

namespace GeometryNamespace {

	/// <summary>
	/// Ctor
	/// A new CS is created on top of the global CS
	/// as ReferenceObject instancee has a strong ownership on the CS instance
	/// and the GlobalCoordSystem is a singleeton class.
	/// GlobalCoordSystem is not owned by any instance
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	ReferenceObject::ReferenceObject(const int theDimensionCount)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		CoordSystem CS = CoordSystem(*GlobalCoordSystem::getGlobalCoordSystem());
		c_referenceCoordSystem = std::shared_ptr<CoordSystem>(&CS);
		inspectDimensionCount(theDimensionCount);
		c_dimensionCount = theDimensionCount;
	}

	/// <summary>
	/// Ctor
	/// Dimension count is DIMENSIONS::D3
	/// </summary>
	ReferenceObject::ReferenceObject(CoordSystem& theReferenceCoordSystem)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
		c_dimensionCount{ DIMENSIONS::D3 },
		c_referenceCoordSystem{ shared_ptr<CoordSystem>(&theReferenceCoordSystem) }
	{
	}

	/// <summary>
	/// The main constructor
	/// </summary>
	ReferenceObject::ReferenceObject(
		const int theDimensionCount,
		CoordSystem& theReferenceCoordSystem)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
		c_referenceCoordSystem{ shared_ptr<CoordSystem>(&theReferenceCoordSystem) }
	{
		inspectDimensionCount(theDimensionCount);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// </summary>
	bool ReferenceObject::operator==(const ReferenceObject& rhs)
	{
		if (&rhs == this) return true;
		if (c_dimensionCount != rhs.getDimensionCount()) return false;
		return *c_referenceCoordSystem == rhs.getReferenceCoordSystem();
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool ReferenceObject::operator!=(const ReferenceObject& rhs) {
		return !operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for ReferenceObject.
	/// Hence, same as == operator
	/// </summary>
	bool ReferenceObject::operator+=(const ReferenceObject& rhs)
	{
		return operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for ReferenceObject.
	/// Hence, same as != operator
	/// </summary>
	bool ReferenceObject::operator-=(const ReferenceObject& rhs)
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	ReferenceObject::~ReferenceObject() {
		Destroy();
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	void ReferenceObject::Destroy() {
		c_dimensionCount = 0;
		c_referenceCoordSystem = nullptr;
	}

	/// <summary>
	/// Used in the copy/move ctor and operators of the child classes
	/// </summary>
	void ReferenceObject::copyBase(const ReferenceObject& rhs)
	{
		GeometryObject::copyBase(rhs);
		c_dimensionCount = rhs.getDimensionCount();
		CoordSystem dummy = rhs.getReferenceCoordSystem();
		c_referenceCoordSystem = shared_ptr<CoordSystem>(&dummy);
	}

	/// <sumrnary>
	/// Returns if the object is 2D
	/// </summary>
	bool ReferenceObject::is2D() const {
		return c_dimensionCount == DIMENSIONS::D2;
	}

	/// <sumrnary>
	/// Static method to determine if the input objects can be members of a 2D object
	/// All inputs shall be 2D and the reference CSs shall be the same
	/// </summary>
	bool ReferenceObject::is2DStrict(
		ARGCOPY(ReferenceObject) theReferenceObject0,
		ARGCOPY(ReferenceObject) theReferenceObject1)
	{
		if (!theReferenceObject0.is2D()) return false;
		if (!theReferenceObject1.is2D()) return false;
		if (theReferenceObject0.getReferenceCoordSystem() == theReferenceObject1.getReferenceCoordSystem()) return true;
		return false;
	}

	/// <sumrnary>
	/// Returns if the object is #D
	/// </summary>
	bool ReferenceObject::is3D() const {
		return c_dimensionCount == DIMENSIONS::D3;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool ReferenceObject::equals(ARGCOPY(ReferenceObject) theReference) const
	{
		return c_referenceCoordSystem->equals(theReference.getReferenceCoordSystem());
	}

	/// <summary>
	/// Equality and geometrical equality are the same for ReferenceObject.
	/// Hence, same as equals(theReferenceObject)
	/// </summary>
	bool ReferenceObject::equalsGeometrically(ARGCOPY(ReferenceObject) theReferenceObject) const
	{
		return equals(theReferenceObject);
	}

	int ReferenceObject::getDimensionCount() const {
		return c_dimensionCount;
	}

	/// <summary>
	/// Getter - Reference CS
	/// </summary>
	CoordSystem& ReferenceObject::getReferenceCoordSystem() const {
		return *c_referenceCoordSystem;
	}

	/// <summary>
	/// Returns the input point if the reference CSs are the same
	/// Otherwise, returns a new point with the same reference CS (global coords are kept the same)
	//! </summary>
	PointBase ReferenceObject::getPointWithMyCoordSystem(ARGCOPY(PointBase) thePoint) const
	{
		if (c_referenceCoordSystem->equals(thePoint.getReferenceCoordSystem()))
		{
			PointBase point { thePoint };
			return point;
		}
		return Point3D(*c_referenceCoordSystem, c_referenceCoordSystem->measurePointCoords(thePoint));
	}

	/// <summary>
	/// Returns the input vector if the reference CSs are the same
	/// Otherwise, returns a new vector with the same reference CS (global components are kept the same)
	/// </summary>
	VectorBase ReferenceObject::getVectorWithMyCoordSystem(ARGCOPY(VectorBase) theVector) const
	{
		if (c_referenceCoordSystem->equals(theVector.getReferenceCoordSystem()))
		{
			return VectorBase(DIMENSIONS::D3, *c_referenceCoordSystem, theVector.getLocalComponents());
		}
		return VectorBase(DIMENSIONS::D3, *c_referenceCoordSystem, c_referenceCoordSystem->measureVectorComponents(theVector));
	}

	/// <summary>
	/// Setter - Dimension count
	/// </summary>
	void ReferenceObject::setDimensionCount(const int theDimensionCount) {
		inspectDimensionCount(theDimensionCount);
		c_dimensionCount = theDimensionCount;
	}

	/// <summary>
	/// Setter - members
	/// Protected methad used in this class and child classes only.
	/// </summary>
	void ReferenceObject::setMembers(const int theDimensionCount, CoordSystem& theReferenceCoordSystem)
	{
		setDimensionCount(theDimensionCount);
		setReferenceCoordSystemBase(theReferenceCoordSystem);
	}

	/// <sumrnary>
	/// Setter - Reference CS (base method)
	/// </summary>
	void ReferenceObject::setReferenceCoordSystemBase(CoordSystem& theReferenceCoordSystem)
	{
		c_referenceCoordSystem = std::shared_ptr<CoordSystem>(&theReferenceCoordSystem);
	}

	/// <sumrnary>
	/// Inspect dimension count
	/// </summary>
	/// <exception> DimensionalityException </exception>
	void ReferenceObject::inspectDimensionCount(const int theDimensionCount) const
	{
		if (theDimensionCount != DIMENSIONS::D2 && theDimensionCount != DIMENSIONS::D3)
		{
			throw DimensionalityException();
		}
	}
}