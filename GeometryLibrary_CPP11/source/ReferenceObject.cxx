// baris.albayrak.ieee@gmail.com

#include "Macros.h"
#include "GeometryObject.hxx"
#include "GeometryMath.hxx"
#include "GeometryParameters.hxx"
#include "GeometryException.hxx"
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
	/// Follows RAII idiom
	/// </summary>
	ReferenceObject::ReferenceObject(const int theDimensionCount)
		:
		GeometryObject(),
		c_dimensionCount { theDimensionCount },
		c_referenceCoordSystem{ GlobalCoordSystem::getGlobalCoordSystem() } { }

	/// <summary>
	/// Ctor
	/// Dimension count is GeometryParameters::DIMENSIONS::D3
	/// </summary>
	ReferenceObject::ReferenceObject(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		:
		GeometryObject(),
		c_referenceCoordSystem{ theReferenceCoordSystem }
	{
	}

	/// <summary>
	/// The main constructor
	/// </summary>
	ReferenceObject::ReferenceObject(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		:
		GeometryObject(),
		c_dimensionCount{ theDimensionCount },
		c_referenceCoordSystem{ theReferenceCoordSystem }
	{
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// </summary>
	bool ReferenceObject::operator==(const ReferenceObject& rhs) const
	{
		if (&rhs == this) return true;
		if (c_dimensionCount != rhs.getDimensionCount()) return false;
		return *c_referenceCoordSystem == *rhs.getReferenceCoordSystem();
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool ReferenceObject::operator!=(const ReferenceObject& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for ReferenceObject.
	/// Hence, same as == operator
	/// </summary>
	bool ReferenceObject::operator+=(const ReferenceObject& rhs) const
	{
		return operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for ReferenceObject.
	/// Hence, same as != operator
	/// </summary>
	bool ReferenceObject::operator-=(const ReferenceObject& rhs) const
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Returns if the object is 2D
	/// </summary>
	bool ReferenceObject::is2D() const {
		return c_dimensionCount == GeometryParameters::DIMENSIONS::D2;
	}

	/// <summary>
	/// Both inputs shall be 2D and the reference CSs shall be the same
	/// </summary>
	bool ReferenceObject::is2DStrict(
		ARGCOPY(ReferenceObject) theReferenceObject0,
		ARGCOPY(ReferenceObject) theReferenceObject1)
	{
		if (!theReferenceObject0.is2D()) return false;
		if (!theReferenceObject1.is2D()) return false;
		if (*theReferenceObject0.getReferenceCoordSystem() == *theReferenceObject1.getReferenceCoordSystem()) return true;
		return false;
	}

	/// <summary>
	/// Returns if the object is #D
	/// </summary>
	bool ReferenceObject::is3D() const {
		return c_dimensionCount == GeometryParameters::DIMENSIONS::D3;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool ReferenceObject::equalsRef(ARGCOPY(ReferenceObject) theReference) const
	{
		return *c_referenceCoordSystem == *theReference.getReferenceCoordSystem();
	}

	/// <summary>
	/// Equality and geometrical equality are the same for ReferenceObject.
	/// Hence, same as equals(theReferenceObject)
	/// </summary>
	bool ReferenceObject::equalsGeometricallyRef(ARGCOPY(ReferenceObject) theReferenceObject) const
	{
		return equalsRef(theReferenceObject);
	}

	int ReferenceObject::getDimensionCount() const {
		return c_dimensionCount;
	}

	/// <summary>
	/// Getter - Reference CS
	/// </summary>
	auto ReferenceObject::getReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>
	{
		return c_referenceCoordSystem;
	}

	/// <summary>
	/// Returns the input point if the reference CSs are the same
	/// Otherwise, returns a new point with the same reference CS (global coords are kept the same)
	//! </summary>
	auto ReferenceObject::getPointWithMyCoordSystem(ARGCOPY(PointBase) thePoint) const -> std::shared_ptr<Point3D>
	{
		if (*c_referenceCoordSystem == *thePoint.getReferenceCoordSystem())
		{
			return std::make_shared<Point3D>(
				c_referenceCoordSystem,
				thePoint.getLocalCoords());
		}
		return std::make_shared<Point3D>(
			c_referenceCoordSystem,
			c_referenceCoordSystem->measurePointCoords(thePoint));
	}

	/// <summary>
	/// Returns the input vector if the reference CSs are the same
	/// Otherwise, returns a new vector with the same reference CS (global components are kept the same)
	/// </summary>
	auto ReferenceObject::getVectorWithMyCoordSystem(ARGCOPY(VectorBase) theVector) const -> std::shared_ptr<Vector3D>
	{
		if (*c_referenceCoordSystem == *theVector.getReferenceCoordSystem())
		{
			return std::make_shared<Vector3D>(
				c_referenceCoordSystem,
				theVector.getLocalComponents());
		}
		return std::make_shared<Vector3D>(
			c_referenceCoordSystem,
			c_referenceCoordSystem->measureVectorComponents(theVector));
	}

	/// <summary>
	/// Setter - Dimension count
	/// </summary>
	void ReferenceObject::setDimensionCount(const int theDimensionCount) {
		c_dimensionCount = theDimensionCount;
	}

	/// <summary>
	/// Setter - members
	/// Protected methad used in this class and child classes only.
	/// </summary>
	void ReferenceObject::setMembersRef(const int theDimensionCount, const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
	{
		setDimensionCount(theDimensionCount);
		setReferenceCoordSystemBase(theReferenceCoordSystem);
	}

	/// <summary>
	/// Setter - Reference CS (base method)
	/// </summary>
	void ReferenceObject::setReferenceCoordSystemBase(const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
	{
		c_referenceCoordSystem = theReferenceCoordSystem;
	}
}
