// baris.albayrak.ieee@gmail.com

#include "GeometryObject.hxx"
#include "ReferenceObject.hxx"
#include "CoordSystem.hxx"

namespace GeometryNamespace {
	long int GeometryObject::c_IDCounter = 0;

	/// <summary>
	/// The default constructor
	/// </summary>
	GeometryObject::GeometryObject() noexcept
	{
		GeometryObject::c_IDCounter++;
		c_ID = c_IDCounter;
	}

	/// <summary>
	/// The constructor with tolerance inputs
	/// </summary>
	GeometryObject::GeometryObject(const double theToleranceGeneral, const double theToleranceSensitive) noexcept
		:
		GeometryObject()
	{
		c_toleranceGeneral = theToleranceGeneral;
		c_toleranceSensitive = theToleranceSensitive;
	}

	/// <summary>
	/// Copy ctor
	/// User defined function is required in order to increment the static ID counter
	/// </summary>
	GeometryObject::GeometryObject(const GeometryObject& rhs)
		:
		c_name { rhs.getName() },
		c_toleranceGeneral { rhs.getToleranceGeneral() },
		c_toleranceSensitive { rhs.getToleranceSensitive() }
	{
		GeometryObject::c_IDCounter++;
		c_ID = GeometryObject::c_IDCounter;
	}

	/// <summary>
	/// Copy assinment
	/// User defined function is required in order to increment the static ID counter
	/// </summary>
	GeometryObject& GeometryObject::operator=(const GeometryObject& rhs) {
		if (&rhs == this) return *this;

		GeometryObject::c_IDCounter++;
		c_ID = GeometryObject::c_IDCounter;
		c_name = rhs.getName();
		c_toleranceGeneral = rhs.getToleranceGeneral();
		c_toleranceSensitive = rhs.getToleranceSensitive();
		return *this;
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
		c_toleranceGeneral = theToleranceGeneral;
	}

	/// <summary>
	/// Setter - Sensitive tolerance
	/// </summary>
	void GeometryObject::setToleranceSensitive(const double& theToleranceSensitive) {
		c_toleranceSensitive = theToleranceSensitive;
	}

	/// <summary>
	/// Inspects the reference CSs of the input ReferenceObject instances.
	/// </summary>
	bool GeometryObject::inspectReferenceCoordSystems(
		ARGCOPY(ReferenceObject) theReference0,
		ARGCOPY(ReferenceObject) theReference1)
	{
		return theReference0.getReferenceCoordSystem() == theReference1.getReferenceCoordSystem();
	}

	/// <summary>
	/// Set the ID with the current static ID (e.g. 42) and create a name with the new ID (e.g. 'Point3D 42')
	/// </summary>
	void GeometryObject::setIDName()
	{
		setID(GeometryObject::c_IDCounter);
		setName(getTypeName(*this));
	}

	/// <summary>
	/// Increment the static ID
	/// </summary>
	void GeometryObject::incrementID()
	{
		GeometryObject::c_IDCounter++;
	}

	/// <summary>
	/// Decrement the static ID
	/// Can be used after creation of a temporary object.
	/// </summary>
	void GeometryObject::decrementID()
	{
		GeometryObject::c_IDCounter--;
	}

	/// <summary>
	/// This method will be used when a new ID and name are required for a newly created object.
	/// This is the only difference comparing with the copy and move ctors and assignments.
	/// Hence, do not use for the temporary objects.
	/// 
	/// A universal reference for the input argument is not used
	/// as it does not make sence to clone a temporary object.
	/// </summary>
	template<typename T>
	T GeometryObject::Clone(const T& arg)
	{
		GeometryObject::c_IDCounter++;
		std::string newName = getTypeName(arg);

		T geometryObject;
		geometryObject.setID(GeometryObject::c_IDCounter);
		geometryObject.setName(newName);
		geometryObject.setToleranceGeneral(arg.getToleranceGeneral());
		geometryObject.setToleranceSensitive(arg.getToleranceSensitive());
		return geometryObject;
	}
}
