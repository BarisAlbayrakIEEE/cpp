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
	IMPLEMENT_STANDARD_RTTIEXT(VectorBase, ReferenceObject)

	/// <summary>
	/// The default constructor for internal usage (protected)
	/// Reference CS is the global CS
	/// Components are all default (i.e. zero magnitude vector, which normally an exception)
	/// </summary>
	VectorBase::VectorBase(const int theDimensionCount)
		:
		ReferenceObject(theDimensionCount),
		c_localComponents{ {} },
		c_magnitude{} { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem) throw (NullptrException)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem)
	{
		c_localComponents = arrayS3{};
		c_magnitude = 0.;
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theLocalComponents) throw (NullptrException)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem)
	{
		inspectLocalComponents(theLocalComponents);
		c_localComponents = theLocalComponents;
		c_magnitude = calculateMagnitude(theLocalComponents);
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theLocalComponents) throw (NullptrException)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem)
	{
		setLocalComponents(theLocalComponents);
	}

	/// <summary>
	/// Ctor
	/// Vector by two points
	/// Reference CS is the reference CS of the points
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	VectorBase::VectorBase(
		const int theDimensionCount,
		ARGCOPY(PointBase) thePoint0,
		ARGCOPY(PointBase) thePoint1) throw (NullptrException, CoordSystemMismatchException)
		: ReferenceObject(theDimensionCount, thePoint0->getReferenceCoordSystem())
	{
		if (!inspectReferenceCoordSystems(thePoint0, thePoint1))
		{
			throw CoordSystemMismatchException(); // CSs
		}

		arrayS3 localComponents = GeometryMath::subtructArrays(thePoint1->getLocalCoords(), thePoint0->getLocalCoords());
		inspectLocalComponents(localComponents);

		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool VectorBase::operator==(const VectorBase& rhs) {
		if (&rhs == this) return true;
		if (ReferenceObject::operator!=(rhs)) return false;
		return GeometryMath::equals(c_localComponents, rhs.getLocalComponents(), c_toleranceGeneral);
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool VectorBase::operator!=(const VectorBase& rhs) {
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool VectorBase::operator+=(const VectorBase& rhs) {
		if (&rhs == this) return true;
		return GeometryMath::equals(getGlobalComponents(), rhs.getGlobalComponents(), c_toleranceGeneral);
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool VectorBase::operator-=(const VectorBase& rhs) {
		return !operator+=(rhs);
	}

	/// <summary>
	/// Vector summation - Member function
	/// </summary>
	OUTVAL(VectorBase) VectorBase::operator+(ARGCOPY(VectorBase) theVector) {
		return add(theVector);
	}

	/// <summary>
	/// Vector subtruction - Member function
	/// </summary>
	OUTVAL(VectorBase) VectorBase::operator-(ARGCOPY(VectorBase) theVector) {
		return subtruct(theVector);
	}

	/// <summary>
	/// Vector summation - Global function
	/// </summary>
	OUTVAL(VectorBase) operator+(ARGCOPY(VectorBase) theVector1, ARGCOPY(VectorBase) theVector2) {
		return theVector1->add(theVector2);
	}

	/// <summary>
	/// Vector subtruction - Global function
	/// </summary>
	OUTVAL(VectorBase) operator-(ARGCOPY(VectorBase) theVector1, ARGCOPY(VectorBase) theVector2) {
		return theVector1->subtruct(theVector2);
	}

	/// <summary>
	/// Use OCCT approach
	/// </summary>
	VectorBase::~VectorBase() {
		Destroy();
	}

	/// <summary>
	/// Use Nullify method of the OCCT Standard_Handle for the object destruction
	/// </summary>
	void VectorBase::Destroy() {
		c_localComponents = arrayS3{ 0., 0., 0. };
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void VectorBase::copyBase(const VectorBase& rhs)
	{
		ReferenceObject::copyBase(rhs);
		std::copy(std::begin(rhs.getLocalComponents()), std::end(rhs.getLocalComponents()), std::begin(c_localComponents));
		c_magnitude = rhs.getMagnitude();
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool VectorBase::equalsBase(ARGCOPY(VectorBase) theVector) const {
		if (theVector.IsNull()) return false;
		return true;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool VectorBase::equals(ARGCOPY(VectorBase) theVector) const {
		if (this == theVector.get()) return true;
		if (!ReferenceObject::equals(theVector)) return false;
		if (!equalsBase(theVector)) return false;
		return GeometryMath::equals(c_localComponents, theVector->getLocalComponents(), c_toleranceGeneral);
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool VectorBase::equalsGeometrically(ARGCOPY(VectorBase) theVector) const
	{
		if (this == theVector.get()) return true;
		if (!equalsBase(theVector)) return false;
		return GeometryMath::equals(getGlobalComponents(), theVector->getGlobalComponents(), c_toleranceGeneral);
	}

	/// <summary>
	/// Inspect local components
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void VectorBase::inspectLocalComponents(const arrayS3& theLocalComponents) const throw (ZeroVectorException)
	{
		if (GeometryMath::equalsZero(theLocalComponents, TOLERANCE_GENERAL))
		{
			throw ZeroVectorException();
		}
	}

	/// <summary>
	/// Getter - Local component - X
	/// </summary>
	double VectorBase::getLocalComponentX() const
	{
		return c_localComponents[0];
	}

	/// <summary>
	/// Getter - Local component - Y
	/// </summary>
	double VectorBase::getLocalComponentY() const
	{
		return c_localComponents[1];
	}

	/// <summary>
	/// Getter - Local component - Z
	/// </summary>
	double VectorBase::getLocalComponentZ() const {
		return c_localComponents[2];
	}

	/// <summary>
	/// Getter - Local components
	/// </summary>
	arrayS3 VectorBase::getLocalComponents() const {
		return c_localComponents;
	}

	/// <summary>
	/// Getter - Global component - X
	/// </summary>
	double VectorBase::getGlobalComponentX() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents[0];

		arrayS3 globalComponents = getGlobalComponents();
		return globalComponents[0];
	}

	/// <summary>
	/// Getter - Global component - Y
	/// </summary>
	double VectorBase::getGlobalComponentY() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents[1];

		arrayS3 globalComponents = getGlobalComponents();
		return globalComponents[1];
	}

	/// <summary>
	/// Getter - Global component - Z
	/// </summary>
	double VectorBase::getGlobalComponentZ() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents[2];

		arrayS3 globalComponents = getGlobalComponents();
		return globalComponents[2];
	}

	/// <summary>
	/// Getter - Global components
	/// </summary>
	arrayS3 VectorBase::getGlobalComponents() const {
		if (c_referenceCoordSystem.IsNull()) return c_localComponents;
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents;

		std::vector<Handle(Vector3D)> axesVectors{ c_referenceCoordSystem->getAxesAsVector() }; // By definition, the reference CS is the global CS
		arrayS3 outGlobalComponents{ arrayS3 {} };
		for (int iAxis = 0; iAxis < DIMENSIONS::D3; iAxis++) {
			arrayS3 componentsVector{};
			try {
				componentsVector = axesVectors[iAxis]->multiply(c_localComponents[iAxis])->getLocalComponents();
			}
			catch (ZeroVectorException) {
				continue;
			}
			outGlobalComponents = GeometryMath::sumArrays(outGlobalComponents, componentsVector);
		}
		return outGlobalComponents;
	}

	/// <summary>
	/// Getter - Slopes
	/// Slope array: Ratio of the components: [y/x, z/y, x/z]
	/// in 2D: [y/x, 0., INF]
	/// </summary>
	arrayS3 VectorBase::getSlopes() const throw (ZeroVectorException)
	{
		// Inspect the local components
		if (
			(GeometryMath::equals(c_localComponents[0], 0., c_toleranceGeneral)) &&
			(GeometryMath::equals(c_localComponents[1], 0., c_toleranceGeneral)) &&
			(GeometryMath::equals(c_localComponents[2], 0., c_toleranceGeneral)))
		{
			throw ZeroVectorException();
		}

		// Determine the slopes
		int componentIndicies[DIMENSIONS::D3][2]{ {0, 1}, {1, 2}, {2, 0} };
		arrayS3 outSlopes{0., 0., GeometryNamespace::INFINITE_VALUE};
		if (is2D()) {
			outSlopes[0] = calculateSlope(c_localComponents[0], c_localComponents[1]);
			return outSlopes;
		}

		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			outSlopes[iCoord] = calculateSlope(componentIndicies[iCoord][0], componentIndicies[iCoord][1]);
		}
		return outSlopes;
	}

	/// <summary>
	/// Getter - Angles
	/// Angle array: Atan of the ratio of the components: [atan(y/x), atan(z/y), atan(x/z)]
	/// in 2D: [atan(y/x), 0., pi / 2]
	/// </summary>
	arrayS3 VectorBase::getAngles() const {
		arrayS3 slopes{ getSlopes() };
		arrayS3 outAngles{ arrayS3{} };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			outAngles[iCoord] = std::atan(slopes[iCoord]);
		}
		return outAngles;
	}

	/// <summary>
	/// Getter - Magnitude
	/// </summary>
	double VectorBase::getMagnitude() const {
		return c_magnitude;
	}

	/// <summary>
	/// Getter - The components of the unit vector in the same direction
	/// </summary>
	arrayS3 VectorBase::getUnitVectorComponents() const {
		return GeometryMath::factorizeArray(c_localComponents, 1. / c_magnitude);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void VectorBase::setMembers(int theDimensionCount, ARGCOPY(CoordSystem) theCoordSystem, const arrayS3& theLocalComponents) {
		ReferenceObject::setMembers(theDimensionCount, theCoordSystem);
		setLocalComponents(theLocalComponents);
	}

	/// <summmary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void VectorBase::setMembers(int theDimensionCount, ARGCOPY(CoordSystem) theCoordSystem, const vectorInput1D& theLocalComponents) {
		ReferenceObject::setMembers(theDimensionCount, theCoordSystem);
		setLocalComponents(theLocalComponents);
	}

	/// <summary>
	/// Setter - Reference CS
	/// Public method for the users.
	/// theKeepGlobalCoordsSame:
	///		true:
	///			Keep the components wrt the GLOBAL CS same and update the components wrt the REFERENCE CS.
	///		false:
	///			Keep the components wrt the REFERENCE CS same and update the components wrt the GLOBAL CS.
	/// </summary>
	void VectorBase::setReferenceCoordSystem(
		ARGCOPY(CoordSystem) theCoordSystem,
		bool theKeepGlobalComponentsSame)
	{
		// No need to update component data
		bool checkCoordSystem{};
		if (!theKeepGlobalComponentsSame) checkCoordSystem = true;
		else if (theCoordSystem->isGlobal()) {
			if (c_referenceCoordSystem->isGlobal()) checkCoordSystem = true;
		}
		else if (theCoordSystem->equals(c_referenceCoordSystem)) checkCoordSystem = true;
		if (checkCoordSystem) {
			ReferenceObject::setReferenceCoordSystemBase(theCoordSystem);
			return;
		}

		// The local components need update
		if (theCoordSystem->isGlobal())
		{
			setLocalComponents(getGlobalComponents());
		}
		else {
			Handle(Point3D) point{ new Point3D(GlobalCoordSystem::getGlobalCoordSystem(), getGlobalComponents()) };
			setLocalComponents(theCoordSystem->measurePointCoords(point));
		}

		// Set the reference CS
		ReferenceObject::setReferenceCoordSystemBase(theCoordSystem);
	}

	/// <summary>
	/// Setter - Local component - X
	/// </summary>
	void VectorBase::setLocalComponentX(const double& theLocalComponentX)
	{
		c_localComponents[0] = theLocalComponentX;
		c_magnitude = calculateMagnitude(c_localComponents);
	}

	/// <summary>
	/// Setter - Local component - Y
	/// </summary>
	void VectorBase::setLocalComponentY(const double& theLocalComponentY)
	{
		c_localComponents[1] = theLocalComponentY;
		c_magnitude = calculateMagnitude(c_localComponents);
	}

	/// <summary>
	/// Setter - Local component - Z
	/// </summary>
	void VectorBase::setLocalComponentZ(const double& theLocalComponentZ)
	{
		c_localComponents[2] = theLocalComponentZ;
		c_magnitude = calculateMagnitude(c_localComponents);
	}

	/// <summary>
	/// Setter - Local components
	/// </summmary>
	void VectorBase::setLocalComponents(const arrayS3& theLocalComponents)
	{
		c_localComponents = theLocalComponents;
		c_magnitude = calculateMagnitude(theLocalComponents);
	}

	/// <summary>
	/// Setter - Local components
	/// </summary>
	void VectorBase::setLocalComponents(const vectorInput1D& theLocalComponents)
	{
		inspectVectorInput(theLocalComponents);
		arrayS3 localComponents = GeometryMath::convertVectorToArray1DS3(theLocalComponents);
		inspectLocalComponents(localComponents);

		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// Calculates the magnitude of the vector
	/// </summary>
	double VectorBase::calculateMagnitude(const arrayS3& theLocalComponents) {
		double magnitude = 0.;
		for (int i = 0; i < DIMENSIONS::D3; i++) {
			magnitude += std::pow(theLocalComponents[i], 2.);
		}
		return std::pow(magnitude, 0.5);
	}

	/// <summary>
	/// Inspects parallelism
	/// Use getToleranceGeneral in order to get the current tolerance value for the calculations
	/// Use setToleranceGeneral before calling this method in order to use another tolerance value
	/// Compares the components (not the slopes) to determine the required condition
	/// Hence, the tolerance for the distance analysis (not sensitive) shall be preferred
	/// Additionally, keep in mind that the parallel vectors may be in reversed directions
	/// Use isInTheSameDirection if the same direction is required
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isParallel(ARGCOPY(VectorBase) theVector) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}
		return isParallel(theVector, c_toleranceGeneral);
	}

	/// <summary>
	/// Inspects parallelism considering the input tolerance
	/// Compares the components (not the slopes) to determine the required condition
	/// Hence, the tolerance for the distance analysis (not sensitive) shall be preferred
	/// Additionally, keep in mind that the parallel vectors may be in reversed directions
	/// Use isInTheSameDirection if the same direction is required
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isParallel(ARGCOPY(VectorBase) theVector, const double& theTolerance) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		arrayS3 globalComponents0{ getGlobalComponents() };
		arrayS3 globalComponents1{ theVector->getGlobalComponents() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (
				!GeometryMath::equals(
					std::fabs(globalComponents0[iCoord]),
					std::fabs(globalComponents1[iCoord]),
					theTolerance)) {
				return false;
			}
		}
		return true;
	}

	/// <summary>
	/// Inspects if in the same direction
	/// Use getToleranceGeneral in order to get the current tolerance va1ue for the calculations
	/// Use setToleranceGeneral before calling this method in order to use another to!erance value
	/// Compares the components (not the slopes) to delerrnine the required condition
	/// Hence, the tolerance for the distance analysis (not sensitive) shall be preferr ed
	/// Additionally, keep in mind that the parallel vectors may be in reversed directions
	/// Use isParallel, if the reversed direction is also acceptable
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isInTheSameDirection(ARGCOPY(VectorBase) theVector) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}
		return isInTheSameDirection(theVector, c_toleranceGeneral);
	}

	/// <summary>
	/// Inspects if in the same direction considering the input tolerance
	/// Compares the components (not the slopes) to determine the required condition
	/// Hence, a general tolerance shall be preferred
	/// Additionally, keep in mind that the parallel vectors may be in reversed directions
	/// Use isParallel, if the reversed direction is also acceptable
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isInTheSameDirection(
		ARGCOPY(VectorBase) theVector,
		const double& theTolerance) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		arrayS3 globalComponents0{ getGlobalComponents() };
		arrayS3 globalComponents1{ theVector->getGlobalComponents() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::equals(globalComponents0[iCoord], globalComponents1[iCoord], theTolerance)) {
				return false;
			}
		}
		return true;
	}

	/// <summary>
	/// Inspects if normal (perpandicular) by measuring the angle between
	/// Use getToleranceSensitive in order to get the current tolerance value for the calculations
	/// Use setToleranceSensitive belore calling this method in order to use another tolerance value
	/// Compares the slopes (not the components) to determine the required condition
	/// Hence, uses the sensitive tolerance value.
	/// Additionally, keep in mind that the parallel vectors may be in reversed directions
	/// Use isParallel, if reversed direction is alsa acceptable
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isNormal(ARGCOPY(VectorBase) theVector) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}
		return isNormal(theVector, c_toleranceSensitive);
	}

	/// <summary>
	/// Inspects if normal (perpandicular) by measuring the angle between considering the input tolerance
	/// Compares the slopes (not the components) to determine the required condition
	/// Hence, a sensitive tolerance shall be preferred
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isNormal(ARGCOPY(VectorBase) theVector, const double& theTolerance) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		double angle{ std::fabs(calculateAngle(theVector)) };
		return bool(
			GeometryMath::equals(angle, M_PI / 2., theTolerance) ||
			GeometryMath::equals(angle, M_PI * 3. / 2., theTolerance));
	}

	/// <summary>
	/// Calculates the angle between the vectors
	/// which by definition inverse of the cosine of the dot product of the vectors
	/// </summary>
	/// <exception> NullptrException </exception>
	double VectorBase::calculateAngle(ARGCOPY(VectorBase) theVector) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		return std::acos(dotProduct(getVectorWithMyCoordSystem(theVector)));
	}

	/// <summary>
	/// Calculates the dot product of the vectors
	/// </summary>
	/// <exception> NullptrException </exception>
	double VectorBase::dotProduct(ARGCOPY(VectorBase) theVector) const throw (NullptrException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		Handle(VectorBase) vector { getVectorWithMyCoordSystem(theVector) };
		arrayS3 localComponents0{ getLocalComponents() };
		arrayS3 localComponents1{ theVector->getLocalComponents() };
		double outDotProduct = 0.;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			outDotProduct += localComponents0[iCoord] * localComponents1[iCoord];
		}
		if (std::fabs(outDotProduct) < c_toleranceGeneral) return 0.;
		return outDotProduct;
	}

	/// <summary>
	/// Returns a 3D vector resulted by the cross product of the vectors
	/// Returns null handle if the two vectors are parallel
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	OUTVAL(Vector3D) VectorBase::crossProduct(ARGCOPY(VectorBase) theVector)
		throw (NullptrException, ZeroVectorException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		arrayS3 localComponents0 = getLocalComponents();
		arrayS3 localComponents1 = getVectorWithMyCoordSystem(theVector)->getLocalComponents();
		arrayS3 vectorComponents{ arrayS3{} };
		vectorComponents[0] = c_localComponents[1] * localComponents1[2] - c_localComponents[2] * localComponents1[1];
		vectorComponents[1] = c_localComponents[2] * localComponents1[0] - c_localComponents[0] * localComponents1[2];
		vectorComponents[2] = c_localComponents[0] * localComponents1[1] - c_localComponents[1] * localComponents1[0];

		Handle(Vector3D) outCrossProduct;
		try {
			outCrossProduct = new Vector3D(vectorComponents);
		}
		catch (...) {
			return Handle(Vector3D)(0);
		}
		return outCrossProduct;
	}

	/// <su m mary>
	/// Returns a 3D point resulted by the translation of the input point by the vector
	/// </summary>
	/// <exception> NullptrException </exception>
	OUTVAL(Point3D) VectorBase::transformPoint(
		ARGCOPY(PointBase) thePoint,
		const double& theFactor) const throw (NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Get the item in my reference CS
		Handle(PointBase) point{ getPointWithMyCoordSystem(thePoint) };

		arrayS3 translation{ GeometryMath::factorizeArray(c_localComponents, theFactor) };
		arrayS3 finalCoords{ GeometryMath::sumArrays(point->getLocalCoords(), translation) };
		return new Point3D(c_referenceCoordSystem, finalCoords);
	}

	/// <summary>
	/// Returns a 3D vector resulted by the sum with the input vector (vector summation)
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	OUTVAL(VectorBase) VectorBase::add(ARGCOPY(VectorBase) theVector) const throw (NullptrException, ZeroVectorException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		// Get the item in my reference CS
		Handle(VectorBase) vector{ getVectorWithMyCoordSystem(theVector) };
		arrayS3 localCoords{ GeometryMath::sumArrays(c_localComponents, vector->getLocalComponents()) };
		try
		{
			if (is2D() && GeometryMath::equals(localCoords[2], 0., c_toleranceGeneral))
			{
				return new VectorBase(DIMENSIONS::D2, c_referenceCoordSystem, localCoords);
			}
			else
			{
				return new VectorBase(DIMENSIONS::D3, c_referenceCoordSystem, localCoords);
			}
		}
		catch (ZeroVectorException)
		{
			throw;
		}
	}

	/// <summary>
	/// Returns a 3D vector resulted by the sum with the reverse of input vector (vector summation)
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	OUTVAL(VectorBase) VectorBase::subtruct(ARGCOPY(VectorBase) theVector) const throw (NullptrException, ZeroVectorException)
	{
		if (theVector.IsNull())
		{
			throw NullptrException();
		}

		Handle(VectorBase) vector{ getVectorWithMyCoordSystem(theVector) };
		arrayS3 localCoords{ GeometryMath::subtructArrays(c_localComponents, vector->getLocalComponents()) };
		try
		{
			if (is2D() && GeometryMath::equals(localCoords[2], 0., c_toleranceGeneral))
			{
				return new VectorBase(DIMENSIONS::D2, c_referenceCoordSystem, localCoords);
			}
			else
			{
				return new VectorBase(DIMENSIONS::D3, c_referenceCoordSystem, localCoords);
			}
		}
		catch (ZeroVectorException)
		{
			throw;
		}
	}

	/// <summary>
	/// Returns a 2D or 3D vector resulted by a scalar multiplication
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	OUTVAL(VectorBase) VectorBase::multiply(const double& theFactor) const throw (ZeroVectorException) {
		if (GeometryMath::equals(theFactor, 0., c_toleranceGeneral)) throw ZeroVectorException();

		arrayS3 finalCoords = GeometryMath::factorizeArray(c_localComponents, theFactor);
		try
		{
			if (is2D() && GeometryMath::equals(finalCoords[2], 0., c_toleranceGeneral))
			{
				return new VectorBase(DIMENSIONS::D2, c_referenceCoordSystem, finalCoords);
			}
			else
			{
				return new VectorBase(DIMENSIONS::D3, c_referenceCoordSystem, finalCoords);
			}
		}
		catch (ZeroVectorException)
		{
			throw;
		}
	}

	/// <summary>
	/// Create unit vector - Global- X
	/// </summary>
	OUTVAL(VectorBase) VectorBase::createUnitVectorX(int theDimensionCount) {
		return new VectorBase(theDimensionCount, GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{ 1., 0., 0. });
	}

	/// <summary>
	/// Create unit vector - Global- Y
	/// </summary>
	OUTVAL(VectorBase) VectorBase::createUnitVectorY(int theDimensionCount)
	{
		return new VectorBase(theDimensionCount, GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{ 0., 1., 0. });
	}

	/// <summary>
	/// Create unit vector - Global- Z
	/// </summary>
	OUTVAL(Vector3D) VectorBase::createUnitVectorZ()
	{
		return new Vector3D(GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{ 0., 0., 1. });
	}

	/// <summary>
	/// Create unit vector - Input CS - X
	/// </summary>
	OUTVAL(VectorBase) VectorBase::createUnitVectorX(int theDimensionCount, ARGCOPY(CoordSystem) theCoordSystem)
	{
		//return new VectorBase(theDimensionCount, theCoordSystem, arrayS3{ 1., 0., 0. });
		if (theDimensionCount == DIMENSIONS::D2)
		{
			return new Vector2D(theCoordSystem, arrayS3{ 1., 0., 0. });
		}
		return new Vector3D(theCoordSystem, arrayS3{ 1., 0., 0. });
	}

	/// <summary>
	/// Create unit vector - Input CS - Y
	/// </summary>
	OUTVAL(VectorBase) VectorBase::createUnitVectorY(int theDimensionCount, ARGCOPY(CoordSystem) theCoordSystem)
	{
		return new VectorBase(theDimensionCount, theCoordSystem, arrayS3{ 0., 1., 0. });
	}

	/// <summary>
	/// Create unit vector - Input CS - Z
	/// </summary>
	OUTVAL(Vector3D) VectorBase::createUnitVectorZ(ARGCOPY(CoordSystem) theCoordSystem)
	{
		return new Vector3D(theCoordSystem, arrayS3{ 0., 0., 1. });
	}

	/// <summary>
	/// Calculates slope for the given two coordinates
	/// slope = coord1 / coord0
	/// slope =
	///     coord0 == 0.
	///         coord1 == 0.
	///             0.
	///         coord1 > 0.
	///             INF
	///         else
	///             -INF
	///     else
	///         coord1 / coord0
	/// </summary>
	double VectorBase::calculateSlope(const double& theCoord0, const double& theCoord1) const {
		if (GeometryMath::equals(theCoord0, 0., c_toleranceGeneral)) {
			if (GeometryMath::equals(theCoord1, 0., c_toleranceGeneral))
			{
				return 0.;
			}
			if (theCoord1 > 0.) {
				return GeometryNamespace::INFINITE_VALUE;
			}
			return -GeometryNamespace::INFINITE_VALUE;
		}
		if (GeometryMath::equals(theCoord1, 0., c_toleranceGeneral)) {
			return 0.;
		}
		return theCoord1 / theCoord0;
	}

	/// <summary>
	/// Calculates angle for the given two coordinates in radians
	/// angle = atan(coord1 / coord0)
	/// angle =
	///     coord0 == 0.
	///         coord1 == 0.
	///             0.
	///         coord1 > 0.
	///             PI / 2.
	///         else
	///             -PI / 2.
	///     else
	///         atan(coord1 / coord0)
	/// </summary>
	double VectorBase::calculateAngle(const double& theCoord0, const double& theCoord1) const {
		return atan(calculateSlope(theCoord0, theCoord1));
	}

	/// <summary>
	/// Find a component that is non-zero based on the input angles.
	/// Find the component that is surely non-zero basedd on the input angles.
	/// Set the component value to 1. for the nonzero component.
	/// Set the other two components based on the input angles.
	/// </summary>
	arrayS3 VectorBase::calculateComponentsFromAngles(const arrayS3& theAngles) const {
		// Inspect the input angles
		arrayS3 outComponents{0., 0., 0.};
		if (
			(GeometryMath::equals(theAngles[0], 0., c_toleranceSensitive)) &&
			(GeometryMath::equals(theAngles[1], 0., c_toleranceSensitive)) &&
			(GeometryMath::equals(theAngles[2], 0., c_toleranceSensitive)))
		{
			return outComponents;
		}

		// Normalize the angles
		arrayS3 normalizedAngles;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			normalizedAngles[iCoord] = GeometryMath::normalizeAngle(theAngles[iCoord]);
		}

		// Calculate slopes
		double slopes[DIMENSIONS::D3];
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			slopes[iCoord] = std::tan(normalizedAngles[iCoord]);
		}

		// Set critical angle values
		double angle90{ M_PI / 2. };
		int count90{ 4 };
		double* criticals90{ new double[count90] { angle90, angle90, angle90 * 3., -angle90 * 3. } };

		// Find a non-zero component using the angles
		bool checkNonZeroComponent{};
		int iCoord{ -1 };
		int nonZeroComponent{};
		while (!checkNonZeroComponent && iCoord < DIMENSIONS::D3 - 1) {
			iCoord++;
			int iAngle{ -1 };
			while (!checkNonZeroComponent && iAngle < count90 - 1) {
				iAngle++;
				if (!GeometryMath::equals(normalizedAngles[iCoord], criticals90[iAngle], c_toleranceSensitive)) {
					checkNonZeroComponent = true;
					nonZeroComponent = iCoord;
				}
			}
		}
		delete[] criticals90;

		// Set the component value to 1. for the non-zero component and the other two to the corresponding slope value
		outComponents[nonZeroComponent] = 1.;
		if (nonZeroComponent == 0) {
			outComponents[1] = slopes[0];
			outComponents[2] = slopes[1] * outComponents[1];
		}
		else if (nonZeroComponent == 1) {
			outComponents[2] = slopes[1];
			outComponents[0] = slopes[2] * outComponents[2];
		}
		else {
			outComponents[0] = slopes[2];
			outComponents[1] = slopes[0] * outComponents[0];
		}
		return outComponents;
	}
}
