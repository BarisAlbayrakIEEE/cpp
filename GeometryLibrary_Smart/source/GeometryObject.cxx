// baris.albayrak.ieee@gmail.com

#include "GeometryObject.hxx"
#include "ReferenceObject.hxx"
#include "CoordSystem.hxx"

namespace GeometryNamespace {

	/// <summary>
	/// The default constructor
	/// </summary>
	GeometryObject::GeometryObject()
		:
		c_name{ "" },
		c_ID{ "" },
		c_toleranceGeneral{ 0. },
		c_toleranceSensitive{ 0. } { }

	/// <summary>
	/// The constructor with tolerance inputs
	/// </summary>
	GeometryObject::GeometryObject(const double theToleranceGeneral, const double theToleranceSensitive)
	{
		inspectTolerance(theToleranceGeneral);
		inspectTolerance(theToleranceSensitive);
		c_toleranceGeneral = theToleranceGeneral;
		c_toleranceSensitive = theToleranceSensitive;
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	GeometryObject::~GeometryObject() {
		c_toleranceGeneral = 0.;
		c_toleranceSensitive = 0.;
		c_name = "";
		c_ID = "";
	}

	/// <summary>
	/// Used in the copy/move ctor and operators of the child classes
	/// </summary>
	void GeometryObject::copyBase(const GeometryObject& rhs)
	{
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
	void GeometryObject::inspectTolerance(const double& theTolerance) const
	{
		if (theTolerance <= 0.) throw ZeroToleranceException();
	}

	/// <sumrnary>
	/// Inspects the reference CSs of the input ReferenceObject instances.
	/// </summary>
	bool GeometryObject::inspectReferenceCoordSystems(
		ARGCOPY(ReferenceObject) theReference0,
		ARGCOPY(ReferenceObject) theReference1)
	{
		return theReference0.getReferenceCoordSystem() == theReference1.getReferenceCoordSystem();
	}

	/// <sumrnary>
	/// Inspects the input std::vector for the size.
	/// The acceptable sizes are:
	///		2 (e.g. local coords of a Point2D object)
	///		3 (e.g. local coords of a Point2D or Point3D objects)
	/// </summary>
	void GeometryObject::inspectVectorInput(const vectorInput1D& theVector) const
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
	void GeometryObject::inspectVectorInputS32(const vectorInput2D& theVector) const
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
