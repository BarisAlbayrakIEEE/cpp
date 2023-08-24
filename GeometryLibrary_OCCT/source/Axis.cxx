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
	IMPLEMENT_STANDARD_RTTIEXT(Axis, GeometryObject)

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Axis::Axis(
		ARGCOPY(PointBase) thePassingPoint,
		ARGCOPY(VectorBase) theDirectionVector)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (thePassingPoint.IsNull())
		{
			throw NullptrException();
		}
		if (theDirectionVector.IsNull())
		{
			throw NullptrException();
		}

		c_passingPoint = thePassingPoint;
		c_directionVector = theDirectionVector;
		updateEC();
	}

	/// <summary>
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	Axis::Axis(
		ARGCOPY(PointBase) thePoint0,
		ARGCOPY(PointBase) thePoint1)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (thePoint0.IsNull())
		{
			throw NullptrException();
		}
		if (thePoint1.IsNull())
		{
			throw NullptrException();
		}

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
		c_passingPoint = thePoint0;
		c_directionVector = vector;
	}

	/// <summary>
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Axis::Axis(
		ARGCOPY(Point2D) thePassingPoint,
		const double& theAngle)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (thePassingPoint.IsNull())
		{
			throw NullptrException();
		}

		Vector2D* vector;
		try {
			vector = new Vector2D(thePassingPoint);
		}
		catch (...) {
			throw;
		}
		c_passingPoint = thePassingPoint;
		c_directionVector = vector;
	}

	/// <summary>
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Axis::Axis(
		ARGCOPY(PointBase) thePassingPoint,
		const arrayS3& theAngles)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (thePassingPoint.IsNull())
		{
			throw NullptrException();
		}

		Vector3D* vector;
		try {
			vector = new Vector3D(thePassingPoint->getReferenceCoordSystem(), theAngles, NULL);
		}
		catch (...) {
			throw;
		}
		c_passingPoint = thePassingPoint;
		c_directionVector = vector;
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
	/// This method inspects final geometrical equality which is actually the coincicience.
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
	/// Use OCCT approach
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
		c_passingPoint = rhs.getPassingPoint();
		c_directionVector = rhs.getDirectionVector();
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
		return ReferenceObject::is2DStrict(c_passingPoint, c_directionVector);
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
	/// Use Nullify method of the OCCT Standard_Handle for the object destruction
	/// </summary>
	void Axis::Destroy()
	{
		c_passingPoint.Nullify();
		c_directionVector.Nullify();
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool Axis::equalsBase(ARGCOPY(Axis) theAxis) const
	{
		if (this == theAxis.get())
		{
			return true;
		}
		if (theAxis.IsNull())
		{
			return false;
		}
		if (!inspectNullEquality(c_passingPoint, theAxis->getPassingPoint()))
		{
			return false;
		}
		if (!inspectNullEquality(c_directionVector, theAxis->getDirectionVector()))
		{
			return false;
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
		if (!c_passingPoint.IsNull() && !c_passingPoint->equals(theAxis->getPassingPoint()))
		{
			return false;
		}
		if (!c_directionVector.IsNull() && !c_directionVector->equals(theAxis->getDirectionVector()))
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
		if (!c_passingPoint.IsNull() && !includes(theAxis->getPassingPoint()))
		{
			return false;
		}
		if (!c_directionVector.IsNull() && !c_directionVector->equalsGeometrically(theAxis->getDirectionVector()))
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
	/// Determines the direction count inspecting if the line has component in z-direction
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
	/// Determines the direction count inspecting if the line has component in z-direction
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
	/// Requires non-null value for the passing point and the direction vector members
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::updateEC()
	{
		if (c_passingPoint.IsNull() || c_directionVector.IsNull())
		{
			throw NullptrException();
		}

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
			c_passingPoint = new Point2D(GlobalCoordSystem::getGlobalCoordSystem(), globalCoords);
			c_directionVector = new Vector2D(GlobalCoordSystem::getGlobalCoordSystem(), globalComponents);
		}
		else {
			c_passingPoint = new Point3D(GlobalCoordSystem::getGlobalCoordSystem(), globalCoords);
			c_directionVector = new Vector3D(GlobalCoordSystem::getGlobalCoordSystem(), globalComponents);
		}
	}

	/// <summary>
	/// Getter - Passing point
	/// </summary>
	OUTVAL(PointBase) Axis::getPassingPoint() const
	{
		return c_passingPoint;
	}

	/// <summary>
	/// Getter - Direction vector
	/// </summary>
	OUTVAL(VectorBase) Axis::getDirectionVector() const
	{
		return c_directionVector;
	}

	/// <summary>
	/// Getter - Equation Coefficients (EC)
	/// </summary>
	arrayS32 Axis::getEC() const
	{
		return c_EC;
	}

	/// <summary>
	/// Creates a point by translating the passing point of the line by the direction vector
	/// </summary>
	OUTVAL(PointBase) Axis::createPoint(const double& theFactor) const
	{
		if (GeometryMath::equals(theFactor, 0., c_toleranceGeneral))
		{
			return c_passingPoint;
		}
		return c_directionVector->transformPoint(c_passingPoint, theFactor);
	}

	/// <summary>
	/// Returns the location vector of the passing point of the line.
	/// Location vector is the vector from the global CS origin to the point
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	OUTVAL(Vector3D) Axis::getPassingPointAsVector() const
	{
		try
		{
			return new Vector3D(
				GlobalCoordSystem::getGlobalCoordSystem(),
				c_passingPoint->getGlobalCoords());
		}
		catch (...)
		{
			throw;
		}
	}

	/// <summary>
	/// Returns the location vector of the point on the line defined as:
	///		Transform the passing point using createPoint method.
	///		Location vector is the vector from the global CS origin to the point
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	OUTVAL(Vector3D) Axis::getPointAsVector(const double& theFactor) const
	{
		Handle(PointBase) point = createPoint(theFactor);
		try
		{
			return new Vector3D(
				GlobalCoordSystem::getGlobalCoordSystem(),
				point->getGlobalCoords());
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
	OUTVAL(CoordSystem) Axis::getCommonReferenceCoordSystem() const
	{
		if (!is2D())
		{
			return Handle(CoordSystem)(0);
		}
		return c_passingPoint->getReferenceCoordSystem();
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
		ARGCOPY(PointBase) thePassingPoint,
		ARGCOPY(VectorBase) theDirectionVector)
	{
		if (thePassingPoint.IsNull())
		{
			throw NullptrException();
		}
		if (theDirectionVector.IsNull())
		{
			throw NullptrException();
		}

		c_passingPoint = thePassingPoint;
		c_directionVector = theDirectionVector;
		updateEC();
	}

	/// <summary>
	/// Setter - Passing point
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setPassingPoint(ARGCOPY(PointBase) thePassingPoint)
	{
		if (thePassingPoint.IsNull())
		{
			throw NullptrException();
		}

		c_passingPoint = thePassingPoint;
		updateEC();
	}

	/// <summary>
	/// Setter - Direction vector
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setDirectionVector(ARGCOPY(VectorBase) theDirectionVector)
	{
		if (theDirectionVector.IsNull())
		{
			throw NullptrException();
		}

		c_directionVector = theDirectionVector;
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
	/// <exception> NullptrException </exception>
	bool Axis::isParallel(ARGCOPY(Axis) theAxis) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return c_directionVector->isParallel(theAxis->getDirectionVector());
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isParallel(
		ARGCOPY(Axis) theAxis,
		const double& theTolerance) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return c_directionVector->isParallel(theAxis->getDirectionVector(), theTolerance);
	}

	/// <summary>
	/// Shortcut to the corresponding melhod of VectorBase
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isInTheSameDirection(ARGCOPY(Axis) theAxis) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return c_directionVector->isInTheSameDirection(theAxis->getDirectionVector());
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isInTheSameDirection(
		ARGCOPY(Axis) theAxis,
		const double& theTolerance) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return c_directionVector->isInTheSameDirection(theAxis->getDirectionVector(), theTolerance);
	}

	/// <summary>
	/// Shortcut to the corresponding method of VectorBase
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isNormal(ARGCOPY(Axis) theAxis) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return c_directionVector->isNormal(theAxis->getDirectionVector());
	}

	/// <summary>
	/// Shortcut o the corresponding method of VectorBase
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isNormal(
		ARGCOPY(Axis) theAxis,
		const double& theTolerance) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return c_directionVector->isNormal(theAxis->getDirectionVector(), theTolerance);
	}

	/// <summary>
	/// Returns if the lines are skew
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isSkew(ARGCOPY(Axis) theAxis)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theAxis) };
		return intersectionResults.first == -1;
	}

	/// <summary>
	/// Returns if the lines are skew
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::isSkew(ARGCOPY(Line) theLine)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theLine->getAxis()) };
		if (
			intersectionResults.first == -1 ||
			(intersectionResults.first == 0 && !theLine->includes(intersectionResults.second))) {
			return true;
		}
		return false;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::includes(ARGCOPY(PointBase) thePoint) const
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Calculate the axis parameter at the point 1ocation
		arrayS3 globalCoords{ thePoint->getGlobalCoords() };
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
	/// <exception> NullptrException </exception>
	bool Axis::intersects(ARGCOPY(Axis) theAxis)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theAxis) };
		return intersectionResults.first == 0;
	}

	/// <summary>
	/// Returns if the line intersects with the input line
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::intersects(ARGCOPY(Line) theLine)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		Handle(Axis) axis1{ theLine->getAxis() };
		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(axis1) };
		if (intersectionResults.first != 0) return false;
		return theLine->includes(intersectionResults.second);
	}

	/// <summary>
	/// Returns if the line coincides with the input line
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::coincides(ARGCOPY(Axis) theAxis)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theAxis) };
		return intersectionResults.first == 1;
	}

	/// <summary>
	/// Returns if the line coincides with the input line
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Axis::coincides(ARGCOPY(Line) theLine)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theLine->getAxis()) };
		return intersectionResults.first == 1;
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input axis
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	/// <exception> NullptrException </exception>
	std::pair<INTERSECTION1, Handle(PointBase)> Axis::intersect(ARGCOPY(Axis) theAxis)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		// Intersection
		bool chkParallelism{};
		Handle(Vector3D) crossProduct = c_directionVector->crossProduct(theAxis->getDirectionVector());
		if (crossProduct.IsNull()) {
			chkParallelism = true;
		}
		if (
			chkParallelism ||
			crossProduct.IsNull() ||
			GeometryMath::equals(crossProduct->getMagnitude(), 0., c_toleranceGeneral))
		{
			if (includes(theAxis->getPassingPoint()))
			{
				return std::pair<INTERSECTION1, Handle(PointBase)>{ INTERSECTION1::Coincides1, Handle(PointBase)(0) };
			}
			return std::pair<INTERSECTION1, Handle(PointBase)>{ INTERSECTION1::Skew1, Handle(PointBase)(0) };
		}
		return intersectBase(theAxis, crossProduct->getGlobalComponents());
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	/// <exception> NullptrException </exception>
	std::pair<INTERSECTION1, Handle(PointBase)> Axis::intersect(ARGCOPY(Line) theLine)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		// The axis
		std::pair<INTERSECTION1, Handle(PointBase)> intersectionResults{ intersect(theLine->getAxis()) };
		if (intersectionResults.first != 0)
			return std::pair<INTERSECTION1, Handle(PointBase)>{ intersectionResults.first, Handle(PointBase)(0) };

		// The line
		intersectionResults.first = (
			theLine->includes(intersectionResults.second)
			? INTERSECTION1::Intersects1
			: INTERSECTION1::Skew1);
		return intersectionResults;
	}

	/// <summary>
	/// Projects the input point onto the line
	/// </summary>
	/// <exception> NullptrException </exception>
	OUTVAL(PointBase) Axis::project(ARGCOPY(PointBase) thePoint) const
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Determine the projection
		Handle(Vector3D) vector0{ new Vector3D(c_passingPoint, thePoint) }; // Vector from the passing point to the point
		double distance{ vector0->dotProduct(c_directionVector) }; // Distance from the passing point to the projection
		Handle(Vector3D) vector1{ new Vector3D(c_passingPoint) }; // The location vector of the passing point
		Handle(VectorBase) vector2{ vector1->add(c_directionVector->multiply(distance)) }; // The location vector of the projection

		if (is2D())
		{
			return new Point2D(vector2->getReferenceCoordSystem(), vector2->getLocalComponents());
		}
		return new Point3D(vector2->getReferenceCoordSystem(), vector2->getLocalComponents());
	}

	/// <summary>
	/// Calculates distance to an axis
	/// </summary>
	/// <exception> NullptrException </exception>
	double Axis::calculateDistance(ARGCOPY(Axis) theAxis)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		// The cross product of the direction vectors
		bool chkParallelism{};
		Handle(VectorBase) crossProduct = c_directionVector->crossProduct(theAxis->getDirectionVector());
		if (crossProduct.IsNull()) {
			chkParallelism = true;
		}

		// The vector between the passing points
		Handle(Vector3D) vectorPassingPoint0{ new Vector3D(c_passingPoint) };
		Handle(Vector3D) vectorPassingPoint1{ new Vector3D(theAxis->getPassingPoint()) };
		Handle(VectorBase) vectorBetweenLines{ vectorPassingPoint0->subtruct(vectorPassingPoint1) };

		// Parallel: d = |V x (R2 - R1)| / |V|	 V, R1 and R2 are vectors and || demonstrates the magnitude
		// where V i s the direction vector, R1 and R 2 are the position vectors of the passing points of the lines
		double outDistance{};
		if (chkParallelism || GeometryMath::equals(crossProduct->getMagnitude(), 0., c_toleranceGeneral))
		{
			return (
				c_directionVector->crossProduct(vectorBetweenLines)->getMagnitude() /
				c_directionVector->getMagnitude());
		}

		// Not parallel: d = (V1 x V2) . (R2 - R1) / |V1 x V2|    V, R1 and R2 are vectors and || demonstrates the magnitude
		// where V is the direction vector, R1 and R2 are the position vectors of the passing points of the lines
		return std::fabs(crossProduct->dotProduct(vectorBetweenLines) / crossProduct->getMagnitude());
	}

	/// <summary>
	/// Calculates distance to a point.
	/// </summary>
	/// <exception> NullptrException </exception>
	double Axis::calculateDistance(ARGCOPY(PointBase) thePoint) const
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		Handle(PointBase) projection{ project(thePoint) };
		return thePoint->calculateDistance(projection);
	}

	/// <summary>
	/// Calculates distance to a line.
	/// Calculate the distance to the axis passing through the line.
	/// The point corresponding to the calculated distance may not be included by the line.
	/// Three cases are possible:
	///		Case 1: The axiss are parallel:
	///			Return the distance to the created axis
	///		Case 2: The axiss intersect:
	///			Return 0. if the intersection of the axiss is included by the line
	///			Otherwise, returns the smaller of the two distances from the intersection to the line end points
	///		Case 3: The axiss are skew:
	///			Finds the closest points on the axiss,
	///			If the line includes the closest point on the 2nd axis, 
	///				returns the distance from the closest point on this line
	///				to the closest point on the 2nd axis
	///			Otherwise, returns the smaller of the two distances from the closest point on this line
	///				to the line end points.
	/// </summary>
	/// <exception> NullptrException </exception>
	double Axis::calculateDistance(ARGCOPY(Line) theLine)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		// The axis
		Handle(Axis) axis1{ theLine->getAxis() };
		double distanceInfinite{ calculateDistance(axis1) };

		// The axiss are parallel
		if (isParallel(axis1)) return distanceInfinite;

		// The axiss intersect
		if (GeometryMath::equals(distanceInfinite, 0., c_toleranceGeneral)) {
			std::pair<int, Handle(PointBase)> intersectionResults{ intersect(axis1) };
			if (theLine->includes(intersectionResults.second)) return 0.;

			std::vector<Handle(PointBase)> endPoints1 { theLine->getEndPoints() };
			double outDistance{ std::fmin(calculateDistance(endPoints1[0]), calculateDistance(endPoints1[1])) };
			return outDistance;
		}

		// The axiss are skew
		std::vector<Handle(PointBase)> closestPoints{ findClosestPoints(axis1) };
		if (theLine->includes(closestPoints[1]))
		{
			double outDistance{ closestPoints[0]->calculateDistance(closestPoints[1]) };
			return outDistance;
		}

		std::vector<Handle(PointBase)> endPoints1 = theLine->getEndPoints();
		double outDistance{ std::fmin(
			closestPoints[0]->calculateDistance(endPoints1[0]),
			closestPoints[0]->calculateDistance(endPoints1[1])) };

		return outDistance;
	}

	/// <summary>
	/// Finds the points on the two axiss which are the closest points.
	/// Solve for t0, t1 and t2 (3 linear equations):
	///	P0 + t0V0 + t2V2 = P1 + t1V1 where capital letters are vectors
	///	V2 = V1 x V0 (cross product)
	///	Ps are the position vectors of the passing points
	///	and Vs are the direction vectors
	/// Returns one of the infiniteLy many point couples in case of parallel lines
	/// </summary>
	/// <exception> NullptrException </exception>
	std::vector<Handle(PointBase)> Axis::findClosestPoints(ARGCOPY(Axis) theAxis)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		// Initialize the output
		std::vector<Handle(PointBase)> outPoints;

		// Get the passing points and the direction vectors
		Handle(PointBase) passingPoint0{ getPassingPoint() };
		Handle(PointBase) passingPoint1{ theAxis->getPassingPoint() };
		Handle(VectorBase) directionVector0{ getDirectionVector() };
		Handle(VectorBase) directionVector1{ theAxis->getDirectionVector() };

		bool chkParallelism{};
		Handle(VectorBase) crossProduct0 = directionVector1->crossProduct(directionVector0);
		if (crossProduct0.IsNull()) {
			chkParallelism = true;
		}

		// Two lines are parallel:
		// Create a line perpandicular to both lines at the passing point of this
		// The passing point and the intersection of the new line with the input line are the outputs
		if (chkParallelism || GeometryMath::equals(crossProduct0->getMagnitude(), 0., c_toleranceGeneral)) {
			outPoints.push_back(passingPoint0);
			Handle(VectorBase) crossProduct1{ crossProduct0->crossProduct(directionVector0)};
			Handle(Axis) axis{ new Axis(passingPoint0, crossProduct1) };
			std::pair<int, Handle(PointBase)> intersectionResults{ intersect(axis)};
			outPoints.push_back(intersectionResults.second);

			return outPoints;
		}

		// Skew lines: P0 + t0V0 + t2V2 = P1 + tM (see the method docstring)
		// Get the coords of the passing points and the components of the direction vectors
		arrayS3 P0{ passingPoint0->getGlobalCoords() };
		arrayS3 P1{ passingPoint1->getGlobalCoords() };
		arrayS3 V0{ directionVector0->getGlobalComponents() };
		arrayS3 V1{ directionVector1->getGlobalComponents() };
		arrayS3 V2{ crossProduct0->getGlobalComponents() };

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
		outPoints[1] = theAxis->createPoint(roots[1]);

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
		Handle(VectorBase) unitVector{ VectorBase::createUnitVectorX(DIMENSIONS::D3) };
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

		Handle(VectorBase) unitVector{ VectorBase::createUnitVectorX(DIMENSIONS::D3) };
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
		Handle(VectorBase) unitVector{ VectorBase::createUnitVectorY(DIMENSIONS::D3) };
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

		Handle(VectorBase) unitVector{ VectorBase::createUnitVectorY(DIMENSIONS::D3) };
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

		Handle(VectorBase) unitVector{ VectorBase::createUnitVectorZ() };
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

		Handle(VectorBase) unitVector{ VectorBase::createUnitVectorZ() };
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
	std::pair<INTERSECTION1, Handle(PointBase)> Axis::intersectBase(
		ARGCOPY(Axis) theAxis,
		const arrayS3& theCrossProductComponents) const
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}

		// Determine the axis parameters at the intersection
		// No need to inspect exception as the parallelism is inspected before
		arrayS32 EC1{ theAxis->getEC() };
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
			return std::pair<INTERSECTION1, Handle(PointBase)>{ INTERSECTION1::Skew1, Handle(PointBase)(0) };
		}

		arrayS3 intersectionCoords;
		intersectionCoords[0] = globalCoord0_x;
		intersectionCoords[1] = globalCoord0_y;
		intersectionCoords[2] = globalCoord0_z;
		Handle(Point3D) intersection = new Point3D(intersectionCoords);

		// Outputs
		Handle(CoordSystem) referenceCoordSystem0 = c_directionVector->getReferenceCoordSystem();
		Handle(CoordSystem) referenceCoordSystem1 = theAxis->getDirectionVector()->getReferenceCoordSystem();
		std::pair<INTERSECTION1, Handle(PointBase)> outIntersectionResults;
		outIntersectionResults.first = INTERSECTION1::Intersects1;
		if (is2D() != theAxis->is2D() || referenceCoordSystem0 != referenceCoordSystem1)
		{
			outIntersectionResults.second = intersection;
		}
		else
		{
			arrayS3 localCoords = referenceCoordSystem0->measurePointCoords(intersection);
			if (is2D() && theAxis->is2D())
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
