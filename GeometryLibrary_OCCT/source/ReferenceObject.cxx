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
	IMPLEMENT_STANDARD_RTTIEXT(ReferenceObject, GeometryObject)

	/// <summary>
	/// Ctor
	/// Reference CS is the global CS
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	ReferenceObject::ReferenceObject(const int theDimensionCount)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
		c_referenceCoordSystem{ GlobalCoordSystem::getGlobalCoordSystem() }
	{
		inspectDimensionCount(theDimensionCount);
		c_dimensionCount = theDimensionCount;
	}

	/// <summary>
	/// Ctor
	/// Dimension count is DIMENSIONS::D3
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	ReferenceObject::ReferenceObject(ARGCOPY(CoordSystem) theReferenceCoordSystem)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
		c_dimensionCount{ DIMENSIONS::D3 }
	{
		if (theReferenceCoordSystem.IsNull())
		{
			throw NullptrException();
		}

		c_referenceCoordSystem = (
			!theReferenceCoordSystem.IsNull() ?
			theReferenceCoordSystem :
			GlobalCoordSystem::getGlobalCoordSystem());
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	ReferenceObject::ReferenceObject(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (theReferenceCoordSystem.IsNull())
		{
			throw NullptrException();
		}
		inspectDimensionCount(theDimensionCount);
		c_dimensionCount = theDimensionCount;
		c_referenceCoordSystem = (
			!theReferenceCoordSystem.IsNull() ?
			theReferenceCoordSystem :
			GlobalCoordSystem::getGlobalCoordSystem());
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// </summary>
	bool ReferenceObject::operator==(const ReferenceObject& rhs)
	{
		if (&rhs == this) return true;
		if (c_dimensionCount != rhs.getDimensionCount()) return false;
		return c_referenceCoordSystem == rhs.getReferenceCoordSystem();
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
	/// Use OCCT approach
	/// </summary>
	ReferenceObject::~ReferenceObject() {
		Destroy();
	}

	/// <summary>
	/// Used in the copy/move ctor and operators of the child classes
	/// </summary>
	void ReferenceObject::copyBase(const ReferenceObject& rhs)
	{
		GeometryObject::copyBase(rhs);
		c_dimensionCount = rhs.getDimensionCount();
		c_referenceCoordSystem = rhs.getReferenceCoordSystem();
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
		if (!inspectNullEquality(theReferenceObject0, theReferenceObject1)) return false;
		if (!theReferenceObject0->is2D()) return false;
		if (!theReferenceObject1->is2D()) return false;

		Handle(CoordSystem) coordSystem1 = theReferenceObject0->getReferenceCoordSystem();
		Handle(CoordSystem) coordSystem2 = theReferenceObject1->getReferenceCoordSystem();
		if (!inspectNullEquality(coordSystem1, coordSystem2)) return false;
		if (coordSystem1 == coordSystem2) return true;
		return false;
	}

	/// <sumrnary>
	/// Returns if the object is #D
	/// </summary>
	bool ReferenceObject::is3D() const {
		return c_dimensionCount == DIMENSIONS::D3;
	}

	/// <summary>
	/// Use Nullify method of the OCCT Standard_Handle for the object destruction
	/// </summary>
	void ReferenceObject::Destroy() {
		c_referenceCoordSystem.Nullify();
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool ReferenceObject::equals(ARGCOPY(ReferenceObject) theReference) const
	{
		if (theReference.IsNull()) return false;
		return c_referenceCoordSystem->equals(theReference->getReferenceCoordSystem());
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
	OUTVAL(CoordSystem) ReferenceObject::getReferenceCoordSystem() const {
		return c_referenceCoordSystem;
	}

	/// <summary>
	/// Returns the input point if the reference CSs are the same
	/// Otherwise, returns a new point with the same reference CS (global coords are kept the same)
	//! </summary>
	OUTVAL(Point3D) ReferenceObject::getPointWithMyCoordSystem(ARGCOPY(PointBase) thePoint) const
		throw (NullptrException)
	{
		if (c_referenceCoordSystem.IsNull()) throw NullptrException();
		if (c_referenceCoordSystem->equals(thePoint->getReferenceCoordSystem()))
		{
			return static_cast<Point3D*>(thePoint.get());
		}
		return new Point3D(c_referenceCoordSystem, c_referenceCoordSystem->measurePointCoords(thePoint));
	}

	/// <summary>
	/// Returns the input vector if the reference CSs are the same
	/// Otherwise, returns a new vector with the same reference CS (global components are kept the same)
	/// </summary>
	OUTVAL(Vector3D) ReferenceObject::getVectorWithMyCoordSystem(ARGCOPY(VectorBase) theVector) const
		throw (NullptrException)
	{
		if (c_referenceCoordSystem.IsNull()) throw NullptrException();
		if (c_referenceCoordSystem->equals(theVector->getReferenceCoordSystem()))
		{
			return static_cast<Vector3D*>(theVector.get());
		}
		return new Vector3D(c_referenceCoordSystem, c_referenceCoordSystem->measureVectorComponents(theVector));
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
	void ReferenceObject::setMembers(const int theDimensionCount, ARGCOPY(CoordSystem) theReferenceCoordSystem)
	{
		setDimensionCount(theDimensionCount);
		setReferenceCoordSystemBase(theReferenceCoordSystem);
	}

	/// <sumrnary>
	/// Setter - Reference CS (base method)
	/// </summary>
	void ReferenceObject::setReferenceCoordSystemBase(ARGCOPY(CoordSystem) theReferenceCoordSystem)
	{
		if (theReferenceCoordSystem.IsNull())
		{
			c_referenceCoordSystem = GlobalCoordSystem::getGlobalCoordSystem();
		}
		c_referenceCoordSystem = theReferenceCoordSystem;
	}

	/// <sumrnary>
	/// Inspect dimension count
	/// </summary>
	/// <exception> DimensionalityException </exception>
	void ReferenceObject::inspectDimensionCount(const int theDimensionCount) const throw (DimensionalityException)
	{
		if (theDimensionCount != DIMENSIONS::D2 && theDimensionCount != DIMENSIONS::D3)
		{
			throw DimensionalityException();
		}
	}
}