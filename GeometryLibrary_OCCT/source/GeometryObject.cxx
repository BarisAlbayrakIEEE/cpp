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
	IMPLEMENT_STANDARD_RTTIEXT(GeometryObject, Standard_Transient)

	/// <summary>
	/// The default constructor
	/// </summary>
	GeometryObject::GeometryObject()
		:
		c_isNull{ true },
		c_name{ "" },
		c_ID{ "" },
		c_toleranceGeneral{ 0. },
		c_toleranceSensitive{ 0. } { }

	/// <summary>
	/// The constructor with tolerance inputs
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	GeometryObject::GeometryObject(const double theToleranceGeneral, const double theToleranceSensitive)
	{
		inspectTolerance(theToleranceGeneral);
		inspectTolerance(theToleranceSensitive);
		c_toleranceGeneral = theToleranceGeneral;
		c_toleranceSensitive = theToleranceSensitive;
	}

	/// <summary>
	/// Use OCCT approach
	/// </summary>
	GeometryObject::~GeometryObject() {
		Destroy();
	}

	/// <summary>
	/// Use Nullify method of the OCCT Standard_Handle for the object destruction
	/// </summary>
	void GeometryObject::Destroy() {
		c_toleranceGeneral = 0.;
		c_toleranceSensitive = 0.;
		c_isNull = true;
		c_name = "";
		c_ID = "";
	}

	/// <summary>
	/// Used in the copy/move ctor and operators of the child classes
	/// </summary>
	void GeometryObject::copyBase(const GeometryObject& rhs)
	{
		c_isNull = rhs.getIsNull();
		c_name = rhs.getName();
		c_ID = rhs.getID();
		c_toleranceGeneral = rhs.getToleranceGeneral();
		c_toleranceSensitive = rhs.getToleranceSensitive();
	}

	/// <summary>
	/// Getter - General tolerance
	/// </summary>
	double GeometryObject::getToleranceGeneral() const {
		return c_toleranceGeneral;
	}

	/// <summary>
	/// Getter - Sensitive tolerance
	/// </summary>
	double GeometryObject::getToleranceSensitive() const {
		return c_toleranceSensitive;
	}

	/// <summary>
	/// Setter - General tolerance
	/// </summary>
	void GeometryObject::setToleranceGeneral(const double& theToleranceGeneral) {
		inspectTolerance(theToleranceGeneral);
		c_toleranceGeneral = theToleranceGeneral;
	}

	/// <summary>
	/// Setter - Sensitive tolerance
	/// </summary>
	void GeometryObject::setToleranceSensitive(const double& theToleranceSensitive) {
		inspectTolerance(theToleranceSensitive);
		c_toleranceSensitive = theToleranceSensitive;
	}

	/// <sumrnary>
	/// Inspects the tolerance input
	/// </summary>
	void GeometryObject::inspectTolerance(const double& theTolerance) const throw (ZeroToleranceException)
	{
		if (theTolerance <= 0.) throw ZeroToleranceException();
	}

	/// <sumrnary>
	/// Inspects the reference CSs of the input ReferenceObject instances.
	/// </summary>
	bool GeometryObject::inspectReferenceCoordSystems(
		ARGCOPY(ReferenceObject) theReference0,
		ARGCOPY(ReferenceObject) theReference1) throw (NullptrException)
	{
		if (theReference0.IsNull() || theReference1.IsNull()) throw NullptrException();
		return theReference0->getReferenceCoordSystem() == theReference1->getReferenceCoordSystem();
	}

	/// <sumrnary>
	/// Inspects if the two instances are:
	///		both null or
	///		both not null
	/// </summary>
	bool GeometryObject::inspectNullEquality(ARGCOPY(GeometryObject) theGeometry0, ARGCOPY(GeometryObject) theGeometry1) {
		if (
			(theGeometry0.IsNull() && !theGeometry1.IsNull()) ||
			(!theGeometry0.IsNull() && theGeometry1.IsNull()))
		{
			return false;
		}
		return true;
	}

	/// <sumrnary>
	/// Inspects the input std::vector for the size.
	/// The acceptable sizes are:
	///		2 (e.g. local coords of a Point2D object)
	///		3 (e.g. local coords of a Point2D or Point3D objects)
	/// </summary>
	void GeometryObject::inspectVectorInput(const vectorInput1D& theVector) const throw (ArraySizeException)
	{
		if (theVector.size() != DIMENSIONS::D2 && theVector.size() != DIMENSIONS::D3)
		{
			throw ArraySizeException();
		}
	}

	/// <sumrnary>
	/// Inspects the two dimensional (i.e. nested) input std::vector
	/// (e.g. equation coefficients (EC) of a line)
	/// </summary>
	void GeometryObject::inspectVectorInputS32(const vectorInput2D& theVector) const throw (ArraySizeException)
	{
		if (theVector.size() != DIMENSIONS::D2 && theVector.size() != DIMENSIONS::D3)
		{
			throw ArraySizeException();
		}
		if (theVector[0].size() != 2)
		{
			throw ArraySizeException();
		}
		if (theVector[1].size() != 2)
		{
			throw ArraySizeException();
		}
		if (theVector.size() == DIMENSIONS::D3 && theVector[2].size() != 2)
		{
			throw ArraySizeException();
		}
	}
}
