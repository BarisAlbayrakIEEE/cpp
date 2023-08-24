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
	/// The main ctor
	/// </summary>
	Axis::Axis(
		PointBase& thePassingPoint,
		VectorBase& theDirectionVector)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		c_directionVector = std::shared_ptr<VectorBase>(&theDirectionVector);
		updateEC();
	}

	/// <summary>
	/// Ctor
	/// </summary>
	Axis::Axis(
		PointBase& thePoint0,
		ARGCOPY(PointBase) thePoint1)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		// Determine if 2D
		int dimensionCount{ DIMENSIONS::D3 };
		if (ReferenceObject::is2DStrict(thePoint0, thePoint1))
		{
			dimensionCount = DIMENSIONS::D2;
		}

		VectorBase* vector;
		try {
			vector = new VectorBase(dimensionCount, thePoint0, thePoint1);
		}
		catch (...) {
			throw;
		}
		c_passingPoint = std::shared_ptr<PointBase>(&thePoint0);
		c_directionVector = std::shared_ptr<VectorBase>(vector);
	}

	/// <summary>
	/// Ctor
	/// </summary>
	Axis::Axis(
		Point2D& thePassingPoint,
		const double& theAngle)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		Vector2D* vector;
		try {
			vector = new Vector2D(thePassingPoint);
		}
		catch (...) {
			throw;
		}
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		c_directionVector = std::shared_ptr<VectorBase>(vector);
	}

	/// <summary>
	/// Ctor
	/// </summary>
	Axis::Axis(
		PointBase& thePassingPoint,
		const arrayS3& theAngles)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		Vector3D* vector;
		try {
			vector = new Vector3D(thePassingPoint.getReferenceCoordSystem(), theAngles, NULL);
		}
		catch (...) {
			throw;
		}
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		c_directionVector = std::shared_ptr<VectorBase>(vector);
	}

	/// <summary>
	/// Ctor using the Equation Coefficients
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Axis::Axis(const arrayS32& theEC)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		try
		{
			setEC(theEC);
		}
		catch (...)
		{
			throw;
		}
	}

	/// <summary>
	/// Ctor using the Equation Coefficients
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Axis::Axis(vectorInput2D theEC)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		try
		{
			setEC(theEC);
		}
		catch (...)
		{
			throw;
		}
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Axis::Axis(const Axis& rhs)
	{
		copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Axis& Axis::operator=(const Axis& rhs)
	{
		if (&rhs == this) return *this;

		copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Axis::Axis(Axis&& rhs) noexcept
	{
		copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Axis& Axis::operator=(Axis&& rhs) noexcept
	{
		if (&rhs == this) return *this;

		copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool Axis::operator==(const Axis& rhs)
	{
		if (&rhs == this) return true;
		return GeometryMath::equals(c_EC, rhs.getEC(), c_toleranceGeneral);
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Axis::operator!=(const Axis& rhs)
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincidence.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Axis::operator+=(const Axis& rhs)
	{
		if (&rhs == this) return true;
		if (!c_directionVector->equalsGeometrically(rhs.getDirectionVector())) return false;
		return includes(rhs.getPassingPoint());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Axis::operator-=(const Axis& rhs)
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	Axis::~Axis()
	{
		Destroy();
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void Axis::copyBase(const Axis& rhs)
	{
		GeometryObject::copyBase(rhs);
		c_EC = GeometryMath::copyNestedArray(rhs.getEC());
		PointBase passingPoint = rhs.getPassingPoint();
		VectorBase directionVector { 
			rhs.getDirectionVector().getDimensionCount(),
			rhs.getDirectionVector().getReferenceCoordSystem(),
			rhs.getDirectionVector().getLocalComponents() };
		c_passingPoint = std::shared_ptr<PointBase>(&passingPoint);
		c_directionVector = std::shared_ptr<VectorBase>(&directionVector);
	}

	/// <summary>
	/// Axis, Line, Circle and Plane types do not have additional 2D and 3D types (e.g. Circle2D)
	/// Hence, these types do not have is2D and is3D methods.
	/// They are assumed 3D by default.
	/// See project main docstring in GeometryObject.hxx for more detailed description.
	/// However, these types can be geometrically 2D.
	/// This method defines the conditions to have a 2D object.
	/// </summary>
	bool Axis::is2D() const
	{
		return ReferenceObject::is2DStrict(*c_passingPoint, *c_directionVector);
	}

	/// <summary>
	/// Returns if the instance is not 2D.
	/// See docstring of is2D for details.
	/// </summary>
	bool Axis::is3D() const
	{
		return !is2D();
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	void Axis::Destroy()
	{
		c_passingPoint = nullptr;
		c_directionVector = nullptr;
	}

	/// <summary>
	/// Obsolute but used. Will be removed in the next revision
	/// </summary>
	bool Axis::equalsBase(ARGCOPY(Axis) theAxis) const
	{
		if (this == &theAxis)
		{
			return true;
		}
		return true;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool Axis::equals(ARGCOPY(Axis) theAxis) const
	{
		if (!equalsBase(theAxis))
		{
			return false;
		}
		if (!c_passingPoint->equals(theAxis.getPassingPoint()))
		{
			return false;
		}
		if (!c_directionVector->equals(theAxis.getDirectionVector()))
		{
			return false;
		}
		return true;
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool Axis::equalsGeometrically(ARGCOPY(Axis) theAxis) const
	{
		if (!equalsBase(theAxis))
		{
			return false;
		}
		if (!includes(theAxis.getPassingPoint()))
		{
			return false;
		}
		if (!c_directionVector->equalsGeometrically(theAxis.getDirectionVector()))
		{
			return false;
		}
		return true;
	}

	/// <sumrnary>
	/// Inspect - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::inspectEC(const arrayS32& theEC)
	{
		if (
			GeometryMath::equals(theEC[0][1], 0., c_toleranceGeneral) &&
			GeometryMath::equals(theEC[1][1], 0., c_toleranceGeneral) &&
			GeometryMath::equals(theEC[2][1], 0., c_toleranceGeneral))
		{
			throw ZeroVectorException();
		}
	}

	/// <sumrnar y>
	/// Inspect - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::inspectEC(const vectorInput2D& theEC)
	{
		inspectVectorInputS32(theEC);
		if (
			GeometryMath::equals(theEC[0][1], 0., c_toleranceGeneral) &&
			GeometryMath::equals(theEC[1][1], 0., c_toleranceGeneral) &&
			GeometryMath::equals(theEC[2][1], 0., c_toleranceGeneral))
		{
			throw ZeroVectorException();
		}
	}

	/// </summary>
	/// Determines the direction count inspecting if the axis has component in z-direction
	/// </summary>
	int Axis::determineDimensionCountUsingEC(const arrayS32& theEC) const
	{
		if (
			GeometryMath::equals(theEC[2][0], 0., c_toleranceGeneral) &&
			GeometryMath::equals(theEC[2][1], 0., c_toleranceGeneral))
		{
			return DIMENSIONS::D2;
		}
		return DIMENSIONS::D3;
	}

	/// </summary>
	/// Determines the direction count inspecting if the axis has component in z-direction
	/// </summary>
	int Axis::determineDimensionCountUsingEC(const vectorInput2D& theEC) const
	{
		if (
			GeometryMath::equals(theEC[2][0], 0., c_toleranceGeneral) &&
			GeometryMath::equals(theEC[2][1], 0., c_toleranceGeneral))
		{
			return DIMENSIONS::D2;
		}
		return DIMENSIONS::D3;
	}

	/// <summary>
	/// Internal method to calculate ECs
	/// </summary>
	void Axis::updateEC()
	{
		arrayS3 globalCoords{ c_passingPoint->getGlobalCoords() };
		arrayS3 globalComponents{ c_directionVector->getGlobalComponents() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			c_EC[iCoord][0] = globalCoords[iCoord];
			c_EC[iCoord][1] = globalComponents[iCoord];
		}
	}

	/// <summary>
	/// Internal method to apply EC (i.e. create passing point and directian vector)
	/// </summary>
	void Axis::applyEC() {
		arrayS3 globalCoords;
		arrayS3 globalComponents;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			globalCoords[iCoord] = c_EC[iCoord][0];
			globalComponents[iCoord] = c_EC[iCoord][1];
		}

		if (
			GeometryMath::equals(c_EC[DIMENSIONS::D3 - 1][0], 0., c_toleranceGeneral) &&
			GeometryMath::equals(c_EC[DIMENSIONS::D3 - 1][1], 0., c_toleranceGeneral))
		{
			globalCoords[DIMENSIONS::D3 - 1] = 0.;
			globalComponents[DIMENSIONS::D3 - 1] = 0.;
			PointBase point = Point2D(*GlobalCoordSystem::getGlobalCoordSystem(), globalCoords);
			VectorBase vector = VectorBase::clone(Vector2D(*GlobalCoordSystem::getGlobalCoordSystem(), globalComponents));
			c_passingPoint = shared_ptr<PointBase>(&point);
			c_directionVector = shared_ptr<VectorBase>(&vector);
		}
		else {
			PointBase point = Point3D(*GlobalCoordSystem::getGlobalCoordSystem(), globalCoords);
			VectorBase vector = VectorBase::clone(Vector3D(*GlobalCoordSystem::getGlobalCoordSystem(), globalComponents));
			c_passingPoint = shared_ptr<PointBase>(&point);
			c_directionVector = shared_ptr<VectorBase>(&vector);
		}
	}

	/// <summary>
	/// Getter - Passing point
	/// </summary>
	PointBase& Axis::getPassingPoint() const
	{
		return *c_passingPoint;
	}

	/// <summary>
	/// Getter - Direction vector
	/// </summary>
	VectorBase& Axis::getDirectionVector() const
	{
		return *c_directionVector;
	}

	/// <summary>
	/// Getter - Equation Coefficients (EC)
	/// </summary>
	arrayS32 Axis::getEC() const
	{
		return c_EC;
	}

	/// <summary>
	/// Creates a point by translating the passing point of the axis by the direction vector
	/// </summary>
	PointBase Axis::createPoint(const double& theFactor) const
	{
		if (GeometryMath::equals(theFactor, 0., c_toleranceGeneral))
		{
			return *c_passingPoint;
		}
		return c_directionVector->transformPoint(*c_passingPoint, theFactor);
	}

	/// <summary>
	/// Returns the location vector of the passing point of the axis.
	/// Location vector is the vector from the global CS origin to the point
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector3D Axis::getPassingPointAsVector() const
	{
		try
		{
			return  Vector3D(*GlobalCoordSystem::getGlobalCoordSystem(), c_passingPoint->getGlobalCoords());
		}
		catch (...)
		{
			throw;
		}
	}

	/// <summary>
	/// Returns the location vector of the point on the axis defined as:
	///		Transform the passing point using createPoint method.
	///		Location vector is the vector from the global CS origin to the point
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Vector3D Axis::getPointAsVector(const double& theFactor) const
	{
		PointBase point = createPoint(theFactor);
		try
		{
			return Vector3D(*GlobalCoordSystem::getGlobalCoordSystem(), point.getGlobalCoords());
		}
		catch (...)
		{
			throw;
		}
	}

	/// <summary>
	/// Returns the reference CS which is common to all ReferenceObject members.
	///		Returns a handle with nullptr if the ReferenceObject members have different reference CSs.
	/// See module docstring in Axis.hxx for the details.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	CoordSystem* Axis::getCommonReferenceCoordSystem() const
	{
		if (!c_passingPoint->getReferenceCoordSystem().equals(c_directionVector->getReferenceCoordSystem()))
		{
			return nullptr;
		}
		return &(c_passingPoint->getReferenceCoordSystem());
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::setMembers(const arrayS32& theEC)
	{
		setEC(theEC);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::setMembers(const vectorInput2D& theEC)
	{
		setEC(theEC);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setMembers(
		PointBase& thePassingPoint,
		VectorBase& theDirectionVector)
	{
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		c_directionVector = std::shared_ptr<VectorBase>(&theDirectionVector);
		updateEC();
	}

	/// <summary>
	/// Setter - Passing point
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setPassingPoint(PointBase& thePassingPoint)
	{
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		updateEC();
	}

	/// <summary>
	/// Setter - Direction vector
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setDirectionVector(VectorBase& theDirectionVector)
	{
		c_directionVector = std::shared_ptr<VectorBase>(&theDirectionVector);
		updateEC();
	}

	/// <summary>
	/// Setter - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::setEC(const arrayS32& theEC)
	{
		try
		{
			inspectEC(theEC);
		}
		catch (...)
		{
			throw;
		}
		c_EC = theEC;
		applyEC();
	}

	/// <summary>
	/// Setter - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::setEC(const vectorInput2D& theEC)
	{
		try
		{
			inspectEC(theEC);
		}
		catch (...)
		{
			throw;
		}
		setEC(GeometryMath::convertVectorToArray2DS32(theEC));
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	bool Axis::isParallel(ARGCOPY(Axis) theAxis) const
	{
		return c_directionVector->isParallel(theAxis.getDirectionVector());
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	bool Axis::isParallel(
		ARGCOPY(Axis) theAxis,
		const double& theTolerance) const
	{
		return c_directionVector->isParallel(theAxis.getDirectionVector(), theTolerance);
	}

	/// <summary>
	/// Shortcut to the corresponding melhod of VectorBase
	/// </summary>
	bool Axis::isInTheSameDirection(ARGCOPY(Axis) theAxis) const
	{
		return c_directionVector->isInTheSameDirection(theAxis.getDirectionVector());
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	bool Axis::isInTheSameDirection(
		ARGCOPY(Axis) theAxis,
		const double& theTolerance) const
	{
		return c_directionVector->isInTheSameDirection(theAxis.getDirectionVector(), theTolerance);
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	bool Axis::isNormal(ARGCOPY(Axis) theAxis) const
	{
		return c_directionVector->isNormal(theAxis.getDirectionVector());
	}

	/// <summary>
	/// Shortcut o the corresponding method of VectorBase
	/// </summary>
	bool Axis::isNormal(
		ARGCOPY(Axis) theAxis,
		const double& theTolerance) const
	{
		return c_directionVector->isNormal(theAxis.getDirectionVector(), theTolerance);
	}

	/// <summary>
	/// Returns if the axes are skew
	/// </summary>
	bool Axis::isSkew(ARGCOPY(Axis) theAxis) const
	{
		std::pair<int, PointBase*> intersectionResults{ intersect(theAxis) };
		return intersectionResults.first == -1;
	}

	/// <summary>
	/// Returns if the axes are skew
	/// </summary>
	bool Axis::isSkew(ARGCOPY(Line) theLine) const
	{
		std::pair<int, PointBase*> intersectionResults{ intersect(theLine.getAxis()) };
		if (
			intersectionResults.first == -1 ||
			(intersectionResults.first == 0 && !theLine.includes(*intersectionResults.second))) {
			return true;
		}
		return false;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Axis::includes(ARGCOPY(PointBase) thePoint) const
	{
		// Calculate the axis parameter at the point 1ocation
		arrayS3 globalCoords{ thePoint.getGlobalCoords() };
		double axisParameter{};
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::equals(c_EC[iCoord][1], 0., c_toleranceGeneral)) {
				axisParameter = (globalCoords[iCoord] - c_EC[iCoord][0]) / c_EC[iCoord][1];
				break;
			}
		}

		// Calculate coordinate and check if they are equal to the point coordinates
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			double globalCoord{ axisParameter * c_EC[iCoord][1] + c_EC[iCoord][0] };
			if (!GeometryMath::equals(globalCoord, globalCoords[iCoord], c_toleranceGeneral)) return false;
		}
		return true;
	}

	/// <summary>
	/// Returns if the axis intersects with the input axis
	/// </summary>
	bool Axis::intersects(ARGCOPY(Axis) theAxis) const
	{
		std::pair<int, PointBase*> intersectionResults{ intersect(theAxis) };
		return intersectionResults.first == 0;
	}

	/// <summary>
	/// Returns if the axis intersects with the input line
	/// </summary>
	bool Axis::intersects(ARGCOPY(Line) theLine) const
	{
		Axis axis1{ theLine.getAxis() };
		std::pair<int, PointBase*> intersectionResults{ intersect(axis1) };
		if (intersectionResults.first != 0) return false;
		return theLine.includes(*intersectionResults.second);
	}

	/// <summary>
	/// Returns if the axis coincides with the input axis
	/// </summary>
	bool Axis::coincides(ARGCOPY(Axis) theAxis) const
	{
		std::pair<int, PointBase*> intersectionResults{ intersect(theAxis) };
		return intersectionResults.first == 1;
	}

	/// <summary>
	/// Returns if the axis coincides with the input line
	/// </summary>
	bool Axis::coincides(ARGCOPY(Line) theLine) const
	{
		std::pair<int, PointBase*> intersectionResults{ intersect(theLine.getAxis()) };
		return intersectionResults.first == 1;
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input axis
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	std::pair<INTERSECTION1, PointBase*> Axis::intersect(ARGCOPY(Axis) theAxis) const
	{
		try {
			Vector3D crossProduct = c_directionVector->crossProduct(theAxis.getDirectionVector());
			return intersectBase(theAxis, crossProduct.getGlobalComponents());
		} catch (ZeroVectorException) {
			if (includes(theAxis.getPassingPoint()))
			{
				return std::pair<INTERSECTION1, PointBase*>{ INTERSECTION1::Coincides1, nullptr };
			}
			return std::pair<INTERSECTION1, PointBase*>{ INTERSECTION1::Skew1, nullptr };
		}
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	std::pair<INTERSECTION1, PointBase*> Axis::intersect(ARGCOPY(Line) theLine) const
	{
		// The axis
		std::pair<INTERSECTION1, PointBase*> intersectionResults{ intersect(theLine.getAxis()) };
		if (intersectionResults.first != 0)
			return std::pair<INTERSECTION1, PointBase*>{ intersectionResults.first, nullptr };

		// The line
		intersectionResults.first = (
			theLine.includes(*intersectionResults.second)
			? INTERSECTION1::Intersects1
			: INTERSECTION1::Skew1);
		return intersectionResults;
	}

	/// <summary>
	/// Projects the input point onto the axis
	/// </summary>
	PointBase Axis::project(ARGCOPY(PointBase) thePoint) const
	{
		// Determine the projection
		Vector3D vector0 { *c_passingPoint, thePoint }; // Vector from the passing point to the point
		double distance { vector0.dotProduct(*c_directionVector) }; // Distance from the passing point to the projection
		Vector3D vector1 { *c_passingPoint }; // The location vector of the passing point
		VectorBase vector2{ vector1.add(c_directionVector->multiply(distance)) }; // The location vector of the projection
		if (is2D())
		{
			return Point2D(vector2.getReferenceCoordSystem(), vector2.getLocalComponents());
		}
		return Point3D(vector2.getReferenceCoordSystem(), vector2.getLocalComponents());
	}

	/// <summary>
	/// Calculates distance to an axis
	/// </summary>
	double Axis::calculateDistance(ARGCOPY(Axis) theAxis) const
	{
		// The vector between the passing points
		Vector3D vectorPassingPoint0{ *c_passingPoint };
		Vector3D vectorPassingPoint1{ theAxis.getPassingPoint() };
		VectorBase vectorBetweenLines = vectorPassingPoint0.subtruct(vectorPassingPoint1);

		try {
			// The cross product of the direction vectors
			VectorBase crossProduct = VectorBase::clone(c_directionVector->crossProduct(theAxis.getDirectionVector()));

			// Not parallel: d = (V1 x V2) . (R2 - R1) / |V1 x V2|    V, R1 and R2 are vectors and || demonstrates the magnitude
			// where V is the direction vector, R1 and R2 are the position vectors of the passing points of the lines
			return std::fabs(crossProduct.dotProduct(vectorBetweenLines) / crossProduct.getMagnitude());
		} catch (ZeroVectorException) {
			return (
				c_directionVector->crossProduct(vectorBetweenLines).getMagnitude() /
				c_directionVector->getMagnitude());

		}
	}

	/// <summary>
	/// Calculates distance to a point.
	/// </summary>
	double Axis::calculateDistance(ARGCOPY(PointBase) thePoint) const
	{
		PointBase projection{ project(thePoint) };
		return thePoint.calculateDistance(projection);
	}

	/// <summary>
	/// Calculates distance to a line.
	/// Calculate the distance to the axis passing through the line.
	/// The point corresponding to the calculated distance may not be included by the line.
	/// Three cases are possible:
	///		Case 1: The axes are parallel:
	///			Return the distance to the created axis
	///		Case 2: The axes intersect:
	///			Return 0. if the intersection of the axes is included by the line
	///			Otherwise, returns the smaller of the two distances from the intersection to the line end points
	///		Case 3: The axes are skew:
	///			Finds the closest points on the axiss,
	///			If the line includes the closest point on the 2nd axis, 
	///				returns the distance from the closest point on this line
	///				to the closest point on the 2nd axis
	///			Otherwise, returns the smaller of the two distances from the closest point on this line
	///				to the line end points.
	/// </summary>
	double Axis::calculateDistance(ARGCOPY(Line) theLine) const
	{
		// The axis
		Axis axis1{ theLine.getAxis() };
		double distanceInfinite{ calculateDistance(axis1) };

		// The axes are parallel
		if (isParallel(axis1)) return distanceInfinite;

		// The axes intersect
		if (GeometryMath::equals(distanceInfinite, 0., c_toleranceGeneral)) {
			std::pair<int, PointBase*> intersectionResults{ intersect(axis1) };
			if (theLine.includes(*intersectionResults.second)) return 0.;

			std::vector<PointBase> endPoints1 { theLine.getEndPoints() };
			double outDistance{ std::fmin(calculateDistance(endPoints1[0]), calculateDistance(endPoints1[1])) };
			return outDistance;
		}

		// The axes are skew
		std::vector<PointBase> closestPoints{ findClosestPoints(axis1) };
		if (theLine.includes(closestPoints[1]))
		{
			double outDistance{ closestPoints[0].calculateDistance(closestPoints[1]) };
			return outDistance;
		}

		std::vector<PointBase> endPoints1 = theLine.getEndPoints();
		double outDistance{ std::fmin(
			closestPoints[0].calculateDistance(endPoints1[0]),
			closestPoints[0].calculateDistance(endPoints1[1])) };

		return outDistance;
	}

	/// <summary>
	/// Finds the points on the two axes which are the closest points.
	/// Solve for t0, t1 and t2 (3 linear equations):
	///	P0 + t0V0 + t2V2 = P1 + t1V1 where capital letters are vectors
	///	V2 = V1 x V0 (cross product)
	///	Ps are the position vectors of the passing points
	///	and Vs are the direction vectors
	/// Returns one of the infiniteLy many point couples in case of parallel lines
	/// </summary>
	std::vector<PointBase> Axis::findClosestPoints(ARGCOPY(Axis) theAxis) const
	{
		// Initialize the output
		std::vector<PointBase> outPoints;

		// Get the passing points and the direction vectors
		PointBase passingPoint0{ getPassingPoint() };
		PointBase passingPoint1{ theAxis.getPassingPoint() };
		VectorBase& directionVector0{ getDirectionVector() };
		VectorBase& directionVector1{ theAxis.getDirectionVector() };
		try {
			VectorBase crossProduct0 = VectorBase::clone(directionVector1.crossProduct(directionVector0));
		} catch (ZeroVectorException) {
			// Two lines are parallel:
			// Create a line perpandicular to both lines at the passing point of this
			// The passing point and the intersection of the new line with the input line are the outputs
			outPoints.push_back(passingPoint0);
			VectorBase dummy { c_directionVector->getDimensionCount(), passingPoint0, passingPoint1 };
			VectorBase crossProduct1 = VectorBase::clone(dummy.crossProduct(directionVector0));
			VectorBase crossProduct2 = VectorBase::clone(directionVector0.crossProduct(crossProduct1));
			Axis axis { passingPoint0, crossProduct2 };
			std::pair<int, PointBase*> intersectionResults{ intersect(axis)};
			outPoints.push_back(*intersectionResults.second);

			return outPoints;
		}

		// Skew lines: P0 + t0V0 + t2V2 = P1 + tM (see the method docstring)
		// Get the coords of the passing points and the components of the direction vectors
		VectorBase crossProduct0 = VectorBase::clone(directionVector1.crossProduct(directionVector0));
		arrayS3 P0{ passingPoint0.getGlobalCoords() };
		arrayS3 P1{ passingPoint1.getGlobalCoords() };
		arrayS3 V0{ directionVector0.getGlobalComponents() };
		arrayS3 V1{ directionVector1.getGlobalComponents() };
		arrayS3 V2{ crossProduct0.getGlobalComponents() };

		// Store the coefficients of the thre equations in a matrix for the solution
		arrayS33 coefficients;
		arrayS3 resultants;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			coefficients[iCoord] = { V0[iCoord], -V1[iCoord], V2[iCoord] };
			resultants[iCoord] = P1[iCoord] - P0[iCoord];
		}

		// Solve the equations
		arrayS33 inverseCoefficients{ GeometryMath::calculateMatrixInverseS33(coefficients, c_toleranceSensitive) };
		arrayS3 roots{ GeometryMath::multiplyMatrixToVectorS33S3(inverseCoefficients, resultants) };
		outPoints[0] = createPoint(roots[0]);
		outPoints[1] = theAxis.createPoint(roots[1]);

		return outPoints;
	}

	/// <summary>
	/// Common code for calculateCoord functions
	/// theIndexCoord0: Index of the axis that the coord value is requested
	/// theIndexCoord1: Index of the axis that the coord value is given
	/// theCoord: The coord value in theIndexCoord1 direction
	/// </summary>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoord(
		int theIndexCoord0,
		int theIndexCoord1,
		const double& theCoord) const
	{
		double p0{ c_EC[theIndexCoord0][0] };
		double p1{ c_EC[theIndexCoord1][0] };
		double V0{ c_EC[theIndexCoord0][1] };
		double V1{ c_EC[theIndexCoord1][1] };
		if (GeometryMath::equals(V1, 0., c_toleranceGeneral)) throw AssymptoticLineException();
		return (theCoord - p1) * V0 / V1 + p0;
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input Y, Requested X
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordX_fromCoordY(const double& theCoordY) const
	{
		VectorBase unitVector{ VectorBase::createUnitVectorX(DIMENSIONS::D3) };
		if (c_directionVector->isParallel(unitVector))
		{
			throw AssymptoticLineException();
		}
		return calculateCoord(0, 1, theCoordY);
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input Z, Requested X
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> DimensionalityException </exception>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordX_fromCoordZ(const double& theCoordZ) const
	{
		if (c_directionVector->is2D())
		{
			throw DimensionalityException();
		}

		VectorBase unitVector{ VectorBase::createUnitVectorX(DIMENSIONS::D3) };
		if (c_directionVector->isParallel(unitVector))
		{
			throw AssymptoticLineException();
		}
		return calculateCoord(0, 2, theCoordZ);
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input X, Requested Y
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordY_fromCoordX(const double& theCoordX) const
	{
		VectorBase unitVector{ VectorBase::createUnitVectorY(DIMENSIONS::D3) };
		if (c_directionVector->isParallel(unitVector))
		{
			throw AssymptoticLineException();
		}
		return calculateCoord(1, 0, theCoordX);
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input Z, Requested Y
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> DimensionalityException </exception>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordY_fromCoordZ(const double& theCoordZ) const
	{
		if (c_directionVector->is2D())
		{
			throw DimensionalityException();
		}

		VectorBase unitVector{ VectorBase::createUnitVectorY(DIMENSIONS::D3) };
		if (c_directionVector->isParallel(unitVector))
		{
			throw AssymptoticLineException();
		}
		return calculateCoord(1, 2, theCoordZ);
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input X, Requested Z
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> DimensionalityException </exception>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordZ_fromCoordX(const double& theCoordX) const
	{
		if (c_directionVector->is2D()) throw DimensionalityException();

		Vector3D unitVector{ VectorBase::createUnitVectorZ() };
		if (c_directionVector->isParallel(unitVector))
		{
			throw AssymptoticLineException();
		}
		return calculateCoord(2, 0, theCoordX);
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input Y, Requested Z
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> DimensionalityException </exception>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordZ_fromCoordY(const double& theCoordY) const
	{
		if (c_directionVector->is2D()) throw DimensionalityException();

		Vector3D unitVector{ VectorBase::createUnitVectorZ() };
		if (c_directionVector->isParallel(unitVector))
		{
			throw AssymptoticLineException();
		}
		return calculateCoord(2, 1, theCoordY);
	}

	/// <summary>
	/// Base method for the intersection with an axis.
	/// Use the parametric equations of the lines
	/// At the intersection, the axis parameter shall be the same for the two lines
	/// Possible values for the output integer:
	///	-1: Lines are skew
	///	0: Lines intersect
	///	1: Lines coincide
	/// CAUTI0N:
	/// This method is not public (e.g. internal method).
	/// Member and the geometric inspections will be performed
	/// in the child class method before calling this method.
	///	Hence, this method does not perform inspection on the input data.
	///	Geometric inspection is the existance of the cross product
	/// </summary>
	/// <exception> NullptrException </exception>
	std::pair<INTERSECTION1, PointBase*> Axis::intersectBase(
		ARGCOPY(Axis) theAxis,
		const arrayS3& theCrossProductComponents) const
	{
		// Determine the axis parameters at the intersection
		// No need to inspect exception as the parallelism is inspected before
		arrayS32 EC1{ theAxis.getEC() };
		double p0_x{ c_EC[0][0] };
		double p0_y{ c_EC[1][0] };
		double p0_z{ c_EC[2][0] };
		double p1_x{ EC1[0][0] };
		double p1_y{ EC1[1][0] };
		double p1_z{ EC1[2][0] };
		double V0_x{ c_EC[0][1] };
		double V0_y{ c_EC[1][1] };
		double V0_z{ c_EC[2][1] };
		double V1_x{ EC1[0][1] };
		double V1_y{ EC1[1][1] };
		double V1_z{ EC1[2][1] };
		double axisParameter0;
		double axisParameter1;
		if (!GeometryMath::equals(theCrossProductComponents[0], 0., c_toleranceGeneral)) {
			axisParameter0 = (V1_z * (p1_y - p0_y) - V1_y * (p1_z - p0_z)) / (V0_y * V1_z - V0_z * V1_y);
			axisParameter1 = (V0_z * (p1_y - p0_y) - V0_y * (p1_z - p0_z)) / (V0_y * V1_z - V0_z * V1_y);
		}
		else if (!GeometryMath::equals(theCrossProductComponents[1], 0., c_toleranceGeneral)) {
			axisParameter0 = (V1_x * (p1_z - p0_z) - V1_x * (p1_x - p0_x)) / (V0_z * V1_x - V0_x * V1_z);
			axisParameter1 = (V0_x * (p1_z - p0_z) - V0_x * (p1_x - p0_x)) / (V0_z * V1_x - V0_x * V1_z);
		}
		else {
			axisParameter0 = (V1_y * (p1_x - p0_x) - V1_x * (p1_y - p0_y)) / (V0_x * V1_y - V0_y * V1_x);
			axisParameter1 = (V0_y * (p1_x - p0_x) - V0_x * (p1_y - p0_y)) / (V0_x * V1_y - V0_y * V1_x);
		}

		// Calculate intersection coordinates
		double globalCoord0_x{ V0_x * axisParameter0 + p0_x };
		double globalCoord0_y{ V0_y * axisParameter0 + p0_y };
		double globalCoord0_z{ V0_z * axisParameter0 + p0_z };
		double globalCoord1_x{ V1_x * axisParameter1 + p1_x };
		double globalCoord1_y{ V1_y * axisParameter1 + p1_y };
		double globalCoord1_z{ V1_z * axisParameter1 + p1_z };
		if (
			!GeometryMath::equals(globalCoord0_x, globalCoord1_x, c_toleranceGeneral) ||
			!GeometryMath::equals(globalCoord0_y, globalCoord1_y, c_toleranceGeneral) ||
			!GeometryMath::equals(globalCoord0_z, globalCoord1_z, c_toleranceGeneral))
		{
			return std::pair<INTERSECTION1, PointBase*>{ INTERSECTION1::Skew1, nullptr };
		}

		arrayS3 intersectionCoords;
		intersectionCoords[0] = globalCoord0_x;
		intersectionCoords[1] = globalCoord0_y;
		intersectionCoords[2] = globalCoord0_z;
		Point3D intersection { intersectionCoords };

		// Outputs
		CoordSystem referenceCoordSystem0 = c_directionVector->getReferenceCoordSystem();
		CoordSystem referenceCoordSystem1 = theAxis.getDirectionVector().getReferenceCoordSystem();
		std::pair<INTERSECTION1, PointBase*> outIntersectionResults;
		outIntersectionResults.first = INTERSECTION1::Intersects1;
		if (is2D() != theAxis.is2D() || referenceCoordSystem0 != referenceCoordSystem1)
		{
			outIntersectionResults.second = &intersection;
		}
		else
		{
			arrayS3 localCoords = referenceCoordSystem0.measurePointCoords(intersection);
			if (is2D() && theAxis.is2D())
			{
				outIntersectionResults.second = new Point2D(
					referenceCoordSystem0,
					localCoords);
			}
			else
			{
				outIntersectionResults.second = new Point3D(
					referenceCoordSystem0,
					localCoords);
			}
		}
		return outIntersectionResults;
	}
}
