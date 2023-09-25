// baris.albayrak.ieee@gmail.com

#include "Macros.h"
#include "GeometryObject.hxx"
#include "GeometryParameters.hxx"
#include "GeometryMath.hxx"
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
	/// The default constructor for internal usage (protected)
	/// Reference CS is the global CS
	/// Components are all default (i.e. zero magnitude vector, which normally an exception)
	/// </summary>
	VectorBase::VectorBase(const int theDimensionCount)
		:
		ReferenceObject(theDimensionCount) {}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem) {	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		const std::array<double, 3>& theLocalComponents)
		: ReferenceObject(theDimensionCount),
		c_localComponents{ theLocalComponents }
	{
		inspectLocalComponents(theLocalComponents);
		c_magnitude = calculateMagnitude(theLocalComponents);
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		const std::vector<double, std::allocator<double>>& theLocalComponents)
		: ReferenceObject(theDimensionCount)
	{
		setLocalComponents(theLocalComponents);
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalComponents)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem),
		c_localComponents{ theLocalComponents }
	{
		inspectLocalComponents(theLocalComponents);
		c_magnitude = calculateMagnitude(theLocalComponents);
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	VectorBase::VectorBase(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalComponents)
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
		ARGCOPY(PointBase) thePoint1)
		: ReferenceObject(theDimensionCount, thePoint0.getReferenceCoordSystem())
	{
		if (!inspectReferenceCoordSystems(thePoint0, thePoint1))
		{
			throw CoordSystemMismatchException(); // CSs
		}

		auto localComponents = GeometryMath::subtructArrays(
			thePoint1.getLocalCoords(),
			thePoint0.getLocalCoords());
		inspectLocalComponents(localComponents);

		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool VectorBase::operator==(const VectorBase& rhs) const
	{
		if (&rhs == this) return true;
		if (ReferenceObject::operator!=(rhs)) return false;
		return GeometryMath::equals(c_localComponents, rhs.getLocalComponents(), getToleranceGeneral());
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool VectorBase::operator!=(const VectorBase& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool VectorBase::operator+=(const VectorBase& rhs) const
	{
		if (&rhs == this) return true;
		return GeometryMath::equals(getGlobalComponents(), rhs.getGlobalComponents(), getToleranceGeneral());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool VectorBase::operator-=(const VectorBase& rhs) const
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool VectorBase::equals(ARGCOPY(VectorBase) theVector) const {
		if (this == &theVector) return true;
		if (!ReferenceObject::equalsRef(theVector)) return false;
		return GeometryMath::equals(c_localComponents, theVector.getLocalComponents(), getToleranceGeneral());
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool VectorBase::equalsGeometrically(ARGCOPY(VectorBase) theVector) const
	{
		if (this == &theVector) return true;
		return GeometryMath::equals(getGlobalComponents(), theVector.getGlobalComponents(), getToleranceGeneral());
	}

	/// <summary>
	/// Inspect local components
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void VectorBase::inspectLocalComponents(const std::array<double, 3>& theLocalComponents) const
	{
		if (GeometryMath::equalsZero(theLocalComponents, getToleranceGeneral()))
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
	auto VectorBase::getLocalComponents() const -> std::array<double, 3>
	{
		return c_localComponents;
	}

	/// <summary>
	/// Getter - Global component - X
	/// </summary>
	double VectorBase::getGlobalComponentX() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents[0];

		auto globalComponents = getGlobalComponents();
		return globalComponents[0];
	}

	/// <summary>
	/// Getter - Global component - Y
	/// </summary>
	double VectorBase::getGlobalComponentY() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents[1];

		auto globalComponents = getGlobalComponents();
		return globalComponents[1];
	}

	/// <summary>
	/// Getter - Global component - Z
	/// </summary>
	double VectorBase::getGlobalComponentZ() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents[2];

		auto globalComponents = getGlobalComponents();
		return globalComponents[2];
	}

	/// <summary>
	/// Getter - Global components
	/// </summary>
	auto VectorBase::getGlobalComponents() const -> std::array<double, 3>
	{
		if (c_referenceCoordSystem->isGlobal()) return c_localComponents;

		auto axesVectors{ c_referenceCoordSystem->getAxesAsVector() }; // By definition, the reference CS is the global CS
		std::array<double, 3> outGlobalComponents = { {} };
		for (int iAxis = 0; iAxis < DIMENSIONS::D3; iAxis++) {
			std::array<double, 3> componentsVector = { {} };
			try {
				componentsVector = axesVectors[iAxis]->multiply(c_localComponents[iAxis])->getLocalComponents();
			}
			catch (ZeroVectorException&) {
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
	auto VectorBase::getSlopes() const -> std::array<double, 3>
	{
		// Inspect the local components
		if (
			(GeometryMath::equals(c_localComponents[0], 0., getToleranceGeneral())) &&
			(GeometryMath::equals(c_localComponents[1], 0., getToleranceGeneral())) &&
			(GeometryMath::equals(c_localComponents[2], 0., getToleranceGeneral())))
		{
			throw ZeroVectorException();
		}

		// Determine the slopes
		std::array<std::array<int, 2>, DIMENSIONS::D3> componentIndicies = { { { {0, 1} }, { { 1, 2} }, {{2, 0}} } };
		std::array<double, 3> outSlopes = { {0., 0., std::numeric_limits<unsigned int>::max()} };
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
	auto VectorBase::getAngles() const {
		auto slopes { getSlopes() };
		std::array<double, 3> outAngles = { {} };
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
	auto VectorBase::getUnitVectorComponents() const -> std::array<double, 3>
	{
		return GeometryMath::factorizeArray(c_localComponents, 1. / c_magnitude);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void VectorBase::setMembers(
		int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theCoordSystem,
		const std::array<double, 3>& theLocalComponents)
	{
		ReferenceObject::setMembersRef(theDimensionCount, theCoordSystem);
		setLocalComponents(theLocalComponents);
	}

	/// <summmary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void VectorBase::setMembers(
		int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theCoordSystem,
		const std::vector<double,
		std::allocator<double>>& theLocalComponents)
	{
		ReferenceObject::setMembersRef(theDimensionCount, theCoordSystem);
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
		const std::shared_ptr<CoordSystem>& theCoordSystem,
		bool theKeepGlobalComponentsSame)
	{
		// No need to update component data
		bool checkCoordSystem{};
		if (!theKeepGlobalComponentsSame) checkCoordSystem = true;
		else if (theCoordSystem->isGlobal()) {
			if (c_referenceCoordSystem->isGlobal()) checkCoordSystem = true;
		}
		else if (theCoordSystem->equals(*c_referenceCoordSystem)) checkCoordSystem = true;
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
			Point3D point{ getGlobalComponents() };
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
	void VectorBase::setLocalComponents(const std::array<double, 3>& theLocalComponents)
	{
		c_localComponents = theLocalComponents;
		c_magnitude = calculateMagnitude(theLocalComponents);
	}

	/// <summary>
	/// Setter - Local components
	/// </summary>
	void VectorBase::setLocalComponents(const std::vector<double, std::allocator<double>>& theLocalComponents)
	{
		std::array<double, 3> localComponents = { {} };
		std::copy(theLocalComponents.begin(), theLocalComponents.end(), localComponents.begin());
		inspectLocalComponents(localComponents);

		c_localComponents = localComponents;
		c_magnitude = calculateMagnitude(localComponents);
	}

	/// <summary>
	/// Calculates the magnitude of the vector
	/// </summary>
	double VectorBase::calculateMagnitude(const std::array<double, 3>& theLocalComponents) const
	{
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
	bool VectorBase::isParallel(ARGCOPY(VectorBase) theVector) const
	{
		return isParallel(theVector, getToleranceGeneral());
	}

	/// <summary>
	/// Inspects parallelism considering the input tolerance
	/// Compares the components (not the slopes) to determine the required condition
	/// Hence, the tolerance for the distance analysis (not sensitive) shall be preferred
	/// Additionally, keep in mind that the parallel vectors may be in reversed directions
	/// Use isInTheSameDirection if the same direction is required
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isParallel(ARGCOPY(VectorBase) theVector, const double& theTolerance) const
	{
		auto globalComponents0{ getGlobalComponents() };
		auto globalComponents1{ theVector.getGlobalComponents() };
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
	bool VectorBase::isInTheSameDirection(ARGCOPY(VectorBase) theVector) const
	{
		return isInTheSameDirection(theVector, getToleranceGeneral());
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
		const double& theTolerance) const
	{
		auto globalComponents0{ getGlobalComponents() };
		auto globalComponents1{ theVector.getGlobalComponents() };
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
	bool VectorBase::isNormal(ARGCOPY(VectorBase) theVector) const
	{
		return isNormal(theVector, getToleranceSensitive());
	}

	/// <summary>
	/// Inspects if normal (perpandicular) by measuring the angle between considering the input tolerance
	/// Compares the slopes (not the components) to determine the required condition
	/// Hence, a sensitive tolerance shall be preferred
	/// </summary>
	/// <exception> NullptrException </exception>
	bool VectorBase::isNormal(ARGCOPY(VectorBase) theVector, const double& theTolerance) const
	{
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
	double VectorBase::calculateAngle(ARGCOPY(VectorBase) theVector) const
	{
		return std::acos(dotProduct(*getVectorWithMyCoordSystem(theVector)));
	}

	/// <summary>
	/// Calculates the dot product of the vectors
	/// </summary>
	/// <exception> NullptrException </exception>
	double VectorBase::dotProduct(ARGCOPY(VectorBase) theVector) const
	{
		auto localComponents0{ getGlobalComponents() };
		auto localComponents1{ theVector.getGlobalComponents() };
		double outDotProduct = 0.;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			outDotProduct += localComponents0[iCoord] * localComponents1[iCoord];
		}
		if (std::fabs(outDotProduct) < getToleranceGeneral()) return 0.;
		return outDotProduct;
	}

	/// <summary>
	/// Returns a 3D vector resulted by the cross product of the vectors
	/// Returns null handle if the two vectors are parallel
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroVectorException </exception>
	auto VectorBase::crossProduct(ARGCOPY(VectorBase) theVector) const -> std::shared_ptr<Vector3D>
	{
		auto localComponents = getVectorWithMyCoordSystem(theVector)->getGlobalComponents();
		std::array<double, 3> vectorComponents = { {} };
		vectorComponents[0] = c_localComponents[1] * localComponents[2] - c_localComponents[2] * localComponents[1];
		vectorComponents[1] = c_localComponents[2] * localComponents[0] - c_localComponents[0] * localComponents[2];
		vectorComponents[2] = c_localComponents[0] * localComponents[1] - c_localComponents[1] * localComponents[0];
		if (GeometryMath::equalsZero<double, 3>(vectorComponents, getToleranceGeneral())) {
			throw ZeroVectorException();
		}

		return std::make_shared<Vector3D>(vectorComponents);
	}

	/// <summary>
	/// Returns a 3D point resulted by the translation of the input point by the vector
	/// </summary>
	/// <exception> NullptrException </exception>
	auto VectorBase::transformPoint(
		ARGCOPY(PointBase) thePoint,
		const double& theFactor) const
		-> std::shared_ptr<Point3D>
	{
		// Get the item in my reference CS
		auto translation{ GeometryMath::factorizeArray(c_localComponents, theFactor) };
		auto finalCoords{
			GeometryMath::sumArrays(
				getPointWithMyCoordSystem(thePoint)->getLocalCoords(),
				translation) };

		return std::make_shared<Point3D>(c_referenceCoordSystem, finalCoords);
	}

	/// <summary>
	/// Create unit vector - Global- Z
	/// </summary>
	auto VectorBase::createUnitVectorZ() -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(std::array<double, 3>{{ 0., 0., 1. }});
	}

	/// <summary>
	/// Create unit vector - Input CS - Z
	/// </summary>
	auto VectorBase::createUnitVectorZ(const std::shared_ptr<CoordSystem>& theCoordSystem) -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(theCoordSystem, std::array<double, 3>{{ 0., 0., 1. }});
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
		if (GeometryMath::equals(theCoord0, 0., getToleranceGeneral())) {
			if (GeometryMath::equals(theCoord1, 0., getToleranceGeneral()))
			{
				return 0.;
			}
			if (theCoord1 > 0.) {
				return std::numeric_limits<unsigned long int>::max();
			}
			return std::numeric_limits<long int>::min();
		}
		if (GeometryMath::equals(theCoord1, 0., getToleranceGeneral())) {
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
	auto VectorBase::calculateComponentsFromAngles(const std::array<double, 3>& theAngles) const -> std::array<double, 3>
	{
		// Inspect the input angles
		std::array<double, 3> outComponents = { {} };
		if (
			(GeometryMath::equals(theAngles[0], 0., getToleranceSensitive())) &&
			(GeometryMath::equals(theAngles[1], 0., getToleranceSensitive())) &&
			(GeometryMath::equals(theAngles[2], 0., getToleranceSensitive())))
		{
			return outComponents;
		}

		// Normalize the angles
		std::array<double, 3> normalizedAngles = { {} };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			normalizedAngles[iCoord] = GeometryMath::normalizeAngle(theAngles[iCoord]);
		}

		// Calculate slopes
		std::array<double, 3> slopes = { {} };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			slopes[iCoord] = std::tan(normalizedAngles[iCoord]);
		}

		// Set critical angle values
		double angle90{ M_PI / 2. };
		std::array<double, 4> criticals90 {{angle90, angle90, angle90 * 3., -angle90 * 3.}};

		// Find a non-zero component using the angles
		bool checkNonZeroComponent{};
		int iCoord{ -1 };
		int nonZeroComponent{};
		while (!checkNonZeroComponent && iCoord < DIMENSIONS::D3 - 1) {
			iCoord++;
			int iAngle{ -1 };
			while (!checkNonZeroComponent && iAngle < 3) {
				iAngle++;
				if (!GeometryMath::equals(normalizedAngles[iCoord], criticals90[iAngle], getToleranceSensitive())) {
					checkNonZeroComponent = true;
					nonZeroComponent = iCoord;
				}
			}
		}

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

	/// <summary>
	/// See GeometryObject::Clone(const T& arg) for descriptions.
	/// </summary>
	template<typename T>
	T VectorBase::Clone(const T& arg)
	{
		T point = GeometryObject::Clone(arg);
		point.setMembers(
			arg.getDimensionCount(),
			arg.getReferenceCoordSystem(),
			arg.getLocalCoords());
		return point;
	}
}
