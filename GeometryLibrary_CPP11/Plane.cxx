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
	/// Ctor using the Equation Coefficients
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Plane::Plane(const std::array<double, 4>& theEC)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(theEC);
	}

	/// <summary>
	/// Ctor using the Equation Coefficients
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Plane::Plane(const std::vector<double, std::allocator<double>>& theEC)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(theEC);
	}

	/// <summary>
	/// Ctor using the passing point and normal vector
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	Plane::Plane(
		const std::shared_ptr<Point3D>& thePassingPoint,
		const std::shared_ptr<Vector3D>& theNormalVector)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(thePassingPoint, theNormalVector);
	}

	/// <summary>
	/// Ctor using the passing point and two in-plane vectors
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Plane::Plane(
		const std::shared_ptr<Point3D>& thePassingPoint,
		ARGCOPY(Vector3D) theInPlaneVector0,
		ARGCOPY(Vector3D) theInPlaneVector1)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
		c_passingPoint{ thePassingPoint },
		c_normalVector{ theInPlaneVector0.crossProduct(theInPlaneVector1) } { }

	/// <summary>
	/// Ctor using three points
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> CoincidenceException </exception>
	/// <exception> ColinearPointsException </exception>
	Plane::Plane(
		const std::shared_ptr<Point3D>& thePassingPoint,
		ARGCOPY(Point3D) thePoint1,
		ARGCOPY(Point3D) thePoint2)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		// Inspect inputs
		if (thePassingPoint->coincides(thePoint1)) throw CoincidenceException();
		if (thePassingPoint->coincides(thePoint2)) throw CoincidenceException();

		Axis axis { thePassingPoint, thePoint1 };
		if (axis.includes(thePoint2)) {
			throw ColinearPointsException();
		}

		// Create two vectors in the plane
		Vector3D vector0{ *thePassingPoint, thePoint1 };
		Vector3D vector1{ *thePassingPoint, thePoint2 };
		c_passingPoint = thePassingPoint;
		c_normalVector = vector0.crossProduct(vector1);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool Plane::operator==(const Plane& rhs) const
	{
		if (&rhs == this) return true;
		if (*c_passingPoint != *rhs.getPassingPoint()) return false;
		return *c_normalVector == *rhs.getNormalVector();
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Plane::operator!=(const Plane& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Plane::operator+=(const Plane& rhs) const
	{
		if (&rhs == this) return true;
		if (!c_normalVector->equalsGeometrically(*rhs.getNormalVector())) return false;
		return includes(*rhs.getPassingPoint());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Plane::operator-=(const Plane& rhs) const
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Axis, Line, Circle and Plane types do not have additional 2D and 3D types (e.g. Circle2D)
	/// Hence, these types do not have is2D and is3D methods.
	/// They are assumed 3D by default.
	/// See project main docstring in GeometryObject.hxx for more detailed description.
	/// However, these types can be geometrically 2D.
	/// This method defines the conditions to have a 2D object.
	/// </summary>
	bool Plane::is2D() const
	{
		if (
			!c_passingPoint->is2D() &&
			!GeometryMath::equals(c_passingPoint->getLocalCoordZ(), 0., getToleranceGeneral()))
		{
			return false;
		}
		if (c_passingPoint->getReferenceCoordSystem()->getAxisAsVectorZ()->isParallel(*c_normalVector))
		{
			return true;
		}
		return false;
	}

	/// <summary>
	/// Returns if the instance is not 2D.
	/// See docstring of is2D for details.
	/// </summary>
	bool Plane::is3D() const
	{
		return !is2D();
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool Plane::equals(ARGCOPY(Plane) thePlane) const
	{
		if (this == &thePlane) return true;
		if (!c_passingPoint->equals(*thePlane.getPassingPoint())) return false;
		if (!c_normalVector->equals(*thePlane.getNormalVector())) return false;
		return true;
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool Plane::equalsGeometrically(ARGCOPY(Plane) thePlane) const
	{
		if (this == &thePlane) return true;

		// The passing points may be different even for eqivalent planes.
		// Hence, unequivalence of passing points does not mean unequivalence of the planes.
		// However, equivalence of passing points and normal vectors means equivalent planes.
		if (c_normalVector->equalsGeometrically(*thePlane.getNormalVector())) return false;
		return includes(*thePlane.getPassingPoint());
	}

	/// <summary>
	/// Inspect Equation Coefficients (EC)
	/// At least one of the 1st three coefficients shall be non-zero.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::inspectEC(const std::array<double, 4>& theEC) const
	{
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::equals(theEC[iCoord], 0., getToleranceGeneral())) {
				return;
			}
		}
		throw ZeroVectorException();
	}

	/// <summary>
	/// Inspect Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::inspectEC(const std::vector<double, std::allocator<double>>& theEC) const
	{
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::equals(theEC[iCoord], 0., getToleranceGeneral())) {
				return;
			}
		}
		throw ZeroVectorException();
	}

	/// </summary>
	/// Getter - Passing point
	/// </summary>
	auto Plane::getPassingPoint() const -> std::shared_ptr<Point3D>
	{
		return c_passingPoint;
	}

	/// <summary>
	/// Getter - Normal vector
	/// </summary>
	auto Plane::getNormalVector() const -> std::shared_ptr<Vector3D>
	{
		return c_normalVector;
	}

	/// <summary>
	/// Getter - EC
	/// </summary>
	auto Plane::getEC() const -> std::array<double, 4>
	{
		std::array<double, 4> EC = { {} };

		double EC3 = 0.;
		auto globalCoords{ c_passingPoint->getGlobalCoords() };
		auto globalComponents{ c_normalVector->getLocalComponents() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			EC[iCoord] = globalComponents[iCoord];
			EC3 -= globalComponents[iCoord] * globalCoords[iCoord];
		}
		EC[3] = EC3;
		return EC;
	}

	/// <summary>
	/// Returns the reference CS which is common to all ReferenceObject members.
	///		Returns a handle with nullptr if the ReferenceObject members have different reference CSs.
	/// See module docstring in Axis.hxx for the details.
	/// </summary>
	auto Plane::getCommonReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>
	{
		if (c_passingPoint->getReferenceCoordSystem() == c_normalVector->getReferenceCoordSystem())
		{
			return c_passingPoint->getReferenceCoordSystem();
		}
		return nullptr;
	}

	/// <summary>
	/// Setter for the members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::setMembers(const std::array<double, 4>& theEC)
	{
		applyEC(theEC);
	}

	/// <summary>
	/// Setter for the members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::setMembers(const std::vector<double, std::allocator<double>>& theEC)
	{
		applyEC(theEC);
	}

	/// <summary>
	/// Setter for the members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void Plane::setMembers(
		const std::shared_ptr<Point3D>& thePassingPoint,
		const std::shared_ptr<Vector3D>& theNormalVector)
	{
		c_passingPoint = thePassingPoint;
		c_normalVector = theNormalVector;
	}

	/// <summary>
	/// Setter - Passing point
	/// Creates a point at the globalorigin if the input point is null.
	/// </summary>
	void Plane::setPassingPoint(const std::shared_ptr<Point3D>& thePassingPoint)
	{
		c_passingPoint = thePassingPoint;
	}

	/// <summary>
	/// Setter - Direction vector
	/// Creates a unit vector (z = 1} at the globalorigin if the input vector is null.
	/// </summary>
	void Plane::setNormalVector(const std::shared_ptr<Vector3D>& theNormalVector)
	{
		c_normalVector = theNormalVector;
	}

	/// <summary>
	/// Returns if the input intersets
	/// </summary>
	bool Plane::intersects(ARGCOPY(Axis) theAxis) const
	{
		return intersect(theAxis).first == 0;
	}

	/// <summary>
	/// Returns if the input intersets
	/// </summary>
	bool Plane::intersects(ARGCOPY(Line) theLine) const
	{
		return intersect(theLine).first == 0 ? true : false;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Plane::includes(ARGCOPY(Point3D) thePoint) const
	{
		// Calculate z-coordinates of the point on the plane having x and y coordinates
		// same as the input point and compare with the point's z-coordinate
		auto globalCoords{ thePoint.getGlobalCoords() };
		double coordX{ globalCoords[0] };
		double coordY{ globalCoords[1] };
		std::pair<bool, double> coordPair{ calculateCoordZ(coordX, coordY) };
		return (
			!coordPair.first ||
			GeometryMath::equals(
				coordPair.second,
				globalCoords[2],
				getToleranceGeneral()));
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Plane::includes(ARGCOPY(Axis) theAxis) const
	{
		// Check if passing point of the line is included
		if (!includes(*theAxis.getPassingPoint())) return false;

		// Check if the plane normal vector and the line direction vector are normal to each other
		return c_normalVector->isNormal(*theAxis.getDirectionVector());
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Plane::includes(ARGCOPY(Line) theLine) const
	{
		// Check if passing point of the line is included
		if (!includes(*theLine.getEndPoint0())) return false;

		// Check if the plane normal vector and the line direction vector are normal to each other
		return c_normalVector->isNormal(*theLine.getDirectionVector());
	}

	/// <summary>
	/// Calculates the coord X
	/// </summary>
	auto Plane::calculateCoordX(const double& theCoordY, const double& theCoordZ) const -> std::pair<bool, double>
	{
		// Initialize outputs
		std::pair<bool, double> outResults{ false, 0. };

		// Get EC
		auto EC{ getEC() };

		// Parameterize inputs
		double coord1{ theCoordY };
		double coord2{ theCoordZ };
		double EC0{ EC[0] };
		double EC1{ EC[1] };
		double EC2{ EC[2] };
		double EC3{ EC[3] };

		// Check if the plane includes the input coordinates
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral())) return outResults;
		outResults.first = true;

		// The output coordinate
		outResults.second = -(EC1 * coord1 + EC2 * coord2 + EC3) / EC0;
		return outResults;
	}

	/// <summary>
	/// Calculates the coord Y
	/// </summary>
	auto Plane::calculateCoordY(const double& theCoordZ, const double& theCoordX) const -> std::pair<bool, double>
	{
		// Initialize outputs
		std::pair<bool, double> outResults{ false, 0. };

		// Get EC
		auto EC{ getEC() };

		// Parameterize inputs
		double coord1{ theCoordZ };
		double coord2{ theCoordX };
		double EC0{ EC[0] };
		double EC1{ EC[1] };
		double EC2{ EC[2] };
		double EC3{ EC[3] };

		// Check if the plane includes the input coordinates
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral())) return outResults;
		outResults.first = true;

		// The output coordinate
		outResults.second = -(EC1 * coord1 + EC2 * coord2 + EC3) / EC0;
		return outResults;
	}

	/// <summary>
	/// Calculates the coord Z
	/// </summary>
	auto Plane::calculateCoordZ(const double& theCoordX, const double& theCoordY) const -> std::pair<bool, double>
	{
		// Initialize outputs
		std::pair<bool, double> outResults{ false, 0. };

		// Get EC
		auto EC{ getEC() };

		// Parameterize inputs
		double coord1{ theCoordX };
		double coord2{ theCoordY };
		double EC0{ EC[0] };
		double EC1{ EC[1] };
		double EC2{ EC[2] };
		double EC3{ EC[3] };

		// Check if the plane includes the input coordinates
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral())) return outResults;
		outResults.first = true;

		// The output coordinate
		outResults.second = -(EC1 * coord1 + EC2 * coord2 + EC3) / EC0;
		return outResults;
	}

	/// <summary>
	/// Calculates the coords X and Y
	/// </summary>
	auto Plane::calculateCoordXY(const double& theCoordZ) const
	{
		// Initialize outputs
		std::pair<bool, std::array<double, 2>> outResults{ false, { {0., 0.} } };

		// Get EC
		auto EC{ getEC() };

		// Parameterize inputs
		double coord2{ theCoordZ };
		double EC0{ EC[0] };
		double EC1{ EC[1] };
		double EC2{ EC[2] };
		double EC3{ EC[3] };

		// Check if the plane includes the input coordinate
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral()) && GeometryMath::equals(EC1, 0., getToleranceGeneral()))
		{
			return outResults;
		}
		outResults.first = true;

		// Check if the plane is dependent only on the st output direction
		if (GeometryMath::equals(EC1, 0., getToleranceGeneral()) && GeometryMath::equals(EC2, 0., getToleranceGeneral()))
		{
			outResults.second[0] = -EC3 / EC0;
			return outResults;
		}

		// Check if the plane is dependent only on the 2nd output direction
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral()) && GeometryMath::equals(EC2, 0., getToleranceGeneral()))
		{
			outResults.second[1] = -EC3 / EC1;
			return outResults;
		}

		// Check if the plane is independent of the st output direction
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral())) {
			outResults.second[1] = -(EC3 + EC2 * coord2) / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 2nd output direction
		if (GeometryMath::equals(EC1, 0., getToleranceGeneral())) {
			outResults.second[0] = -(EC3 + EC2 * coord2) / EC0;
			return outResults;
		}

		// All exceptional cases inspected.
		// All directions are dependent
		// Keep the 1st output direction with zero initial value and calculate the 2nd direction correspondingly
		outResults.second[1] = -(EC3 + EC2 * coord2) / EC1;
		return outResults;
	}

	/// <summary>
	/// Calculates the coords Y and Z
	/// </summary>
	auto Plane::calculateCoordYZ(const double& theCoordX) const
	{
		// Initialize outputs
		std::pair<bool, std::array<double, 2>> outResults{ false, { {0., 0.} } };

		// Get EC
		auto EC{ getEC() };

		// Parameteri ze inputs
		double coord2{ theCoordX };
		double EC0{ EC[1] };
		double EC1{ EC[2] };
		double EC2{ EC[0] };
		double EC3{ EC[3] };

		// Check if the plane includes the input coordinate
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral()) && GeometryMath::equals(EC1, 0., getToleranceGeneral()))
		{
			return outResults;
		}
		outResults.first = true;

		// Check if the plane is dependent only on the st output direction
		if (GeometryMath::equals(EC1, 0., getToleranceGeneral()) && GeometryMath::equals(EC2, 0., getToleranceGeneral()))
		{
			outResults.second[0] = -EC3 / EC0;
			return outResults;
		}

		// Check if the plane is dependent only on the 2nd output direction
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral()) && GeometryMath::equals(EC2, 0., getToleranceGeneral()))
		{
			outResults.second[1] = -EC3 / EC1;
			return outResults;
		}

		// Check if the plane is independent of the st output direction
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral())) {
			outResults.second[1] = -(EC3 + EC2 * coord2) / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 2nd output direction
		if (GeometryMath::equals(EC1, 0., getToleranceGeneral())) {
			outResults.second[0] = -(EC3 + EC2 * coord2) / EC0;
			return outResults;
		}

		// All exceptional cases inspected.
		// All directions are dependent
		// Keep the st output direction with zero initial value and calculatethe 2nd direction correspondingly
		outResults.second[1] = -(EC3 + EC2 * coord2) / EC1;
		return outResults;
	}

	/// <summary>
	/// Calculates the coords Z and X
	/// </summary>
	auto Plane::calculateCoordZX(const double& theCoordY) const
	{
		// Initialize outputs
		std::pair<bool, std::array<double, 2>> outResults{ false, { {0., 0.} } };

		// Get EC
		auto EC{ getEC() };

		// Parameterize inputs
		double coord2{ theCoordY };
		double EC0{ EC[2] };
		double EC1{ EC[0] };
		double EC2{ EC[1] };
		double EC3{ EC[3] };

		// Check if the plane includes the input coordinate
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral()) && GeometryMath::equals(EC1, 0., getToleranceGeneral()))
		{
			return outResults;
		}
		outResults.first = true;

		// Check if the plane is dependent only on the st output direction
		if (GeometryMath::equals(EC1, 0., getToleranceGeneral()) && GeometryMath::equals(EC2, 0., getToleranceGeneral()))
		{
			outResults.second[0] = -EC3 / EC0;
			return outResults;
		}

		// Check if the plane is dependent only on the 2nd output direction
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral()) && GeometryMath::equals(EC2, 0., getToleranceGeneral()))
		{
			outResults.second[1] = -EC3 / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 1st output direction
		if (GeometryMath::equals(EC0, 0., getToleranceGeneral())) {
			outResults.second[1] = -(EC3 + EC2 * coord2) / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 2nd output direction
		if (GeometryMath::equals(EC1, 0., getToleranceGeneral())) {
			outResults.second[0] = -(EC3 + EC2 * coord2) / EC0;
			return outResults;
		}

		// All exceptional cases inspected.
		// All directions are dependent
		// Keep the 1st output direction with zero initial value and calculatethe 2nd direction correspondingly
		outResults.second[1] = -(EC3 + EC2 * coord2) / EC1;
		return outResults;
	}

	/// <summary>
	/// Calculates distance to the input plane
	/// </summary>
	double Plane::calculateDistance(ARGCOPY(Plane) thePlane) const
	{
		if (!c_normalVector->isParallel(*thePlane.c_normalVector)) return 0.;
		return calculateDistance(*thePlane.c_passingPoint);
	}

	/// <summary>
	/// Calculates distance to the input point
	/// </summary>
	double Plane::calculateDistance(ARGCOPY(Point3D) thePoint) const
	{
		// Get EC
		auto EC{ getEC() };

		Point3D point { thePoint.getReferenceCoordSystem(), thePoint.getLocalCoords() };
		double value0{};
		double value1{};
		auto globalCoords{ point.getGlobalCoords() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			value0 += EC[iCoord] * globalCoords[iCoord];
			value1 += std::pow(EC[iCoord], 2.);
		}

		value0 += EC[3];
		value1 = std::pow(value1, 0.5);
		return std::fabs(value0 / value1);
	}

	/// <summary>
	/// Calculates distance to the input axis
	/// </summary>
	double Plane::calculateDistance(ARGCOPY(Axis) theAxis) const
	{
		// Check if the two intersects
		if (!c_normalVector->isNormal(*theAxis.getDirectionVector())) return 0.;

		// Find distance from the line's passing point to the plane
		return calculateDistance(*theAxis.getPassingPoint());
	}

	/// <summary>
	/// Calculates distance to the input line
	/// </summary>
	double Plane::calculateDistance(ARGCOPY(Line) theLine) const
	{
		// Check if the two intersects
		if (intersects(theLine)) return 0.;

		// Find distances for the two points
		return std::fmin(
			calculateDistance(*theLine.getEndPoint0()), 
			calculateDistance(*theLine.getEndPoint1()));
	}

	/// <summary>
	/// Returns the intersection status and the intersection object (if exists) with the input plane
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Include
	/// 
	/// Approach:
	/// 1. Determine the direction of the intersection line by cross product of the planes' normal.
	/// 2. Determine the passing point by solving the plane equations together.
	///		The problem has 3 unknowns (x, y and z) but two equations.
	///	Hence, we need to cancel one of the unknowns by assigning a trivial value.
	///	However, this should be done carefully as the intersection line may be independent of one of x, y and z.
	///	This can be determined by inspecting the planes' ECs.
	/// </summary>
	auto Plane::intersect(ARGCOPY(Plane) thePlane) const -> std::pair<INTERSECTION2, std::shared_ptr<Axis>>
	{
		// Check parallelism
		if (c_normalVector->isParallel(*thePlane.getNormalVector()))
		{
			if (GeometryMath::equals(calculateDistance(*thePlane.c_passingPoint), 0., getToleranceGeneral()))
			{
				return std::pair<INTERSECTION2, std::shared_ptr<Axis>>{ INTERSECTION2::Includes2, nullptr };
			}
			return std::pair<INTERSECTION2, std::shared_ptr<Axis>>{ INTERSECTION2::Skew2, nullptr };
		}

		// Determine a point on the intersection line
		// intersection is a line
		// Hence, it passes through all coordinates (including zero) for at least one direction
		// Method:
		//	1. Assume x-coordinate is zero
		//	2. Find y-coordinate by equating the z-coordinate in two plane equations
		//	3. Find z-ccordinate by using one of the plane equations with the known x and y coordinates
		//		if the solution in Step 2 is achieved
		//	4. Otherwise, start from Step 1 by setting y-coordinate to zero and repeat all steps
		//	5. Repeat all steps with z = 0 if Step 4 also fails.
		//        Failure for all three directions cannot happen as the planes are not parallel.
		double coordX{};
		double coordY{};
		double coordZ{};
		double divisor{};
		double denominator{};

		// Get EC
		auto EC1{ getEC() };
		auto EC2{ thePlane.getEC() };

		// Trial 1: x = o.
		coordX = 0.;
		denominator = EC1[1] * EC2[2] - EC2[1] * EC1[2];
		if (!GeometryMath::equals(denominator, 0., getToleranceGeneral())) {
			divisor = EC1[2] * EC2[3] - EC2[2] * EC1[3];
			coordY = divisor / denominator;
		}
		else {
			// Trial 2: y = 0.
			coordY = 0.;
			denominator = EC1[2] * EC2[0] - EC2[2] * EC1[0];
			if (!GeometryMath::equals(denominator, 0., getToleranceGeneral())) {
				divisor = EC1[0] * EC2[3] - EC2[0] * EC1[3];
				coordZ = divisor / denominator;
			}
			else {
				// Trial 3: z = 0.
				coordZ = 0.;
				denominator = EC1[0] * EC2[1] - EC2[0] * EC1[1];
				if (!GeometryMath::equals(denominator, 0., getToleranceGeneral())) {
					divisor = EC1[1] * EC2[3] - EC2[1] * EC1[3];
					coordX = divisor / denominator;
				}
			}
		}
		Point3D passingPoint{ std::array<double, 3>{{ coordX, coordY, coordZ } } };

		// The direction vector of the intersection line
		// No need to inspect exception as the parallelism is checl<ed at the beginning of the method
		auto directionVector{ c_normalVector->crossProduct(*thePlane.c_normalVector) };

		// intersection line
		return std::pair<INTERSECTION2, std::shared_ptr<Axis>>{
			INTERSECTION2::Intersects2,
			std::make_shared<Axis>(std::make_shared<Point3D>(passingPoint), directionVector) };
	}

	/// <summary>
	/// Returns the intersection status and the intersection object (if exists) with the input axis
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// 
	/// Approach:
	/// See wikipedia: Line-plane intersection
	/// </summary>
	auto Plane::intersect(ARGCOPY(Axis) theAxis) const -> std::pair<INTERSECTION2, std::shared_ptr<Point3D>>
	{
		// Check parallelism
		double angle{ std::fabs(c_normalVector->calculateAngle(*theAxis.getDirectionVector())) };
		if (
			GeometryMath::equals(angle, M_PI / 2., getToleranceSensitive()) ||
			GeometryMath::equals(angle, M_PI * 3. / 2., getToleranceSensitive()))
		{
			if (GeometryMath::equals(calculateDistance(*theAxis.getPassingPoint()), 0., getToleranceGeneral()))
			{
				return std::pair<INTERSECTION2, std::shared_ptr<Point3D>>{ INTERSECTION2::Includes2, nullptr };
			}
			return std::pair<INTERSECTION2, std::shared_ptr<Point3D>>{ INTERSECTION2::Skew2, nullptr };
		}

		// Solve the plane's scalar equation (Ax + By + Cz + D = 0) and line's parametric equation
		// (x - x0/Vx = y - y0/Vy = z - z0/Vz = t) together to obtain line's parameter (t)
		// A * (x0 + Vx * t) + 8 * (y0 + Vy * t) + C * (z0 + Vz * t) + D = 0 where t is the only unknown
		// No need to inspect exception as the parallelism is inspected before
		auto EC1{ getEC() };
		auto EC2{ theAxis.getEC() };
		double axisParameter{ -(
			(
				EC1[0] * EC2[0][0] +
				EC1[1] * EC2[0][1] +
				EC1[2] * EC2[0][2] +
				EC1[3]) /
			(
				EC1[0] * EC2[1][0] +
				EC1[0] * EC2[1][1] +
				EC1[0] * EC2[1][2])) };

		// Determine point coords using line1 s parameteric equation with the parameter
		std::array<double, 3> pointCoords{{
			EC2[0][0] + EC2[1][0] * axisParameter,
			EC2[0][1] + EC2[1][1] * axisParameter,
			EC2[0][2] + EC2[1][2] * axisParameter }};

		// Create point
		return std::pair<INTERSECTION2, std::shared_ptr<Point3D>>{
			INTERSECTION2::Intersects2,
				std::make_shared<Point3D>(pointCoords) };
	}

	/// <summary>
	/// Returns the intersection status and the intersection object (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	auto Plane::intersect(ARGCOPY(Line) theLine) const -> std::pair<INTERSECTION2, std::shared_ptr<Point3D>>
	{
		// Inspect intersection for the infinite tine corresponding to the line
		auto intersectionResults{ intersect(*theLine.getAxis()) };
		if (intersectionResults.first != 0)
		{
			return intersectionResults;
		}

		// Check if the intersection point is included by the line
		intersectionResults.first = (
			theLine.includes(*intersectionResults.second)
			? INTERSECTION2::Intersects2
			: INTERSECTION2::Skew2);
		return intersectionResults;
	}

	/// <summary>
	/// Projects the input point onto the plane
	/// </summary>
	auto Plane::project(ARGCOPY(Point3D) thePoint) const -> std::shared_ptr<Point3D>
	{
		// Check if the point is on the plane
		if (includes(thePoint)) {
			return std::make_shared<Point3D>(
				thePoint.getReferenceCoordSystem(),
				thePoint.getLocalCoords());
		}

		// Create a line passing through the input point and in the direction parallel to the planes normal vector
		auto point{ thePoint };
		Axis axis { std::shared_ptr<Point3D>(&point), c_normalVector };

		// The projection is the intersection point
		auto intersectionResults{ intersect(axis) };
		return intersectionResults.second;
	}

	/// <summary>
	/// Projects the input vector onto the plane
	/// </summary>
	auto Plane::project(ARGCOPY(Vector3D) theVector) const -> std::shared_ptr<Vector3D>
	{
		// Create an arbitray axis using the vector
		Point3D arbitraryPoint{ std::array<double, 3>{{}} };
		auto transformPoint { theVector.transformPoint(arbitraryPoint, 2.) };
		Axis axis { std::shared_ptr<Point3D>(&arbitraryPoint), *transformPoint };

		// Project the axis onto the plane
		auto projection { project(axis) };
		if (!projection) return nullptr;

		// The direction vector of the projection is the projected vector
		return projection->getDirectionVector();
	}

	/// <summary>
	/// Projects the input axis onto the plane
	/// </summary>
	auto Plane::project(ARGCOPY(Axis) theAxis) const -> std::shared_ptr<Axis>
	{
		// Take two points on the line and project onto the plane
		auto pointOnPlane0 { *project(*theAxis.getPassingPoint()) };
		auto pointOnLine0 { theAxis.createPoint(2.) };
		Point3D pointOnLine1 { pointOnLine0->getReferenceCoordSystem(), pointOnLine0->getLocalCoords() };
		auto pointOnPlane1 { project(pointOnLine1) };
		if (pointOnPlane0.coincides(*pointOnPlane1))
		{
			return nullptr;
		}
		return std::make_shared<Axis>(std::shared_ptr<Point3D>(&pointOnPlane0), *pointOnPlane1);
	}

	/// <summary>
	/// Projects the input line onto the plane
	/// </summary>
	auto Plane::project(ARGCOPY(Line) theLine) const -> std::shared_ptr<Line>
	{
		// Take two points on the line and project onto the plane
		auto pointOnPlane0 { project(*theLine.getEndPoint0()) };
		auto pointOnPlane1 { project(*theLine.getEndPoint1()) };
		if (pointOnPlane0->coincides(*pointOnPlane1))
		{
			return nullptr;
		}

		// Create line using the projection points
		return std::make_shared<Line>(pointOnPlane0, pointOnPlane1);
	}

	/// <summary>
	/// Creates a point using Equation Coefficients (EC)
	/// </summary>
	/// <exception> UncaughtException </exception>
	auto Plane::createPoint(const double& theFactor) const
	{
		auto EC = getEC();

		// Set a value to one of the thre coordinates and determine the other two coordinates
		double coordX{};
		bool exists{};
		double coordY{};
		double coordZ{};
		auto globalCoords{ c_passingPoint->getGlobalCoords() };
		if (!GeometryMath::equals(EC[0], 0., getToleranceGeneral())) {
			coordX = globalCoords[0] + theFactor;
			auto coordPair{ calculateCoordYZ(coordX) };
			exists = coordPair.first;
		}
		else if (!GeometryMath::equals(EC[1], 0., getToleranceGeneral())) {
			coordY = globalCoords[1] + theFactor;
			auto coordPair{ calculateCoordZX(coordY) };
			exists = coordPair.first;
		}
		else if (!GeometryMath::equals(EC[2], 0., getToleranceGeneral())) {
			coordZ = globalCoords[2] + theFactor;
			auto coordPair{ calculateCoordXY(coordZ) };
			exists = coordPair.first;
		}
		else {
			coordX = globalCoords[0] + theFactor;
			auto coordPair{ calculateCoordYZ(coordX) };
			exists = coordPair.first;
		}
		if (!exists) throw UncaughtException();

		return std::make_shared<Point3D>(std::array<double, 3>{{ coordX, coordY, coordZ}});
	}

	/// <summary>
	/// Apply EC (i.e. create passing point and normal vector)
	/// </summary>
	void Plane::applyEC(const std::array<double, 4>& theEC)
	{
		inspectEC(theEC);

		// Passing point
		double coordX{};
		double coordY{};
		double coordZ{};
		if (std::fabs(theEC[1]) <= getToleranceGeneral() && std::fabs(theEC[2]) <= getToleranceGeneral())
		{
			coordX = -theEC[3] / theEC[0];
		}
		else if (std::fabs(theEC[0]) <= getToleranceGeneral() && std::fabs(theEC[2]) <= getToleranceGeneral())
		{
			coordY = -theEC[3] / theEC[1];
		}
		else if (std::fabs(theEC[0]) <= getToleranceGeneral() && std::fabs(theEC[1]) <= getToleranceGeneral())
		{
			coordZ = -theEC[3] / theEC[2];
		}
		else if (std::fabs(theEC[0]) <= getToleranceGeneral())
		{
			coordZ = -theEC[3] / theEC[2]; // for y = 0
		}
		else if (std::fabs(theEC[1]) <= getToleranceGeneral())
		{
			coordX = -theEC[3] / theEC[0]; // for z = 0
		}
		else if (std::fabs(theEC[2]) <= getToleranceGeneral())
		{
			coordY = -theEC[3] / theEC[1]; // for x = 0
		}
		else
		{
			coordZ = -theEC[3] / theEC[2]; // for x = 0 and y = 0
		}
		Point3D passingPoint{ std::array<double, 3>{{ coordX, coordY, coordZ } } };
		c_passingPoint = std::shared_ptr<Point3D>(&passingPoint);

		// The normal vector
		std::array<double, 3> globalComponents = { {} };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			globalComponents[iCoord] = theEC[iCoord];
		}
		Vector3D normalVector{ globalComponents };
		c_normalVector = std::shared_ptr<Vector3D>(&normalVector);
	}

	/// <summary>
	/// Apply EC (i.e. create passing point and normal vector)
	/// </summary>
	void Plane::applyEC(const std::vector<double, std::allocator<double>>& theEC)
	{
		inspectEC(theEC);

		// Passing point
		double coordX{};
		double coordY{};
		double coordZ{};
		if (std::fabs(theEC[1]) <= getToleranceGeneral() && std::fabs(theEC[2]) <= getToleranceGeneral())
		{
			coordX = -theEC[3] / theEC[0];
		}
		else if (std::fabs(theEC[0]) <= getToleranceGeneral() && std::fabs(theEC[2]) <= getToleranceGeneral())
		{
			coordY = -theEC[3] / theEC[1];
		}
		else if (std::fabs(theEC[0]) <= getToleranceGeneral() && std::fabs(theEC[1]) <= getToleranceGeneral())
		{
			coordZ = -theEC[3] / theEC[2];
		}
		else if (std::fabs(theEC[0]) <= getToleranceGeneral())
		{
			coordZ = -theEC[3] / theEC[2]; // for y = 0
		}
		else if (std::fabs(theEC[1]) <= getToleranceGeneral())
		{
			coordX = -theEC[3] / theEC[0]; // for z = 0
		}
		else if (std::fabs(theEC[2]) <= getToleranceGeneral())
		{
			coordY = -theEC[3] / theEC[1]; // for x = 0
		}
		else
		{
			coordZ = -theEC[3] / theEC[2]; // for x = 0 and y = 0
		}
		Point3D passingPoint{ std::array<double, 3>{{ coordX, coordY, coordZ } } };
		c_passingPoint = std::shared_ptr<Point3D>(&passingPoint);

		// The normal vector
		std::array<double, 3> globalComponents = { {} };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			globalComponents[iCoord] = theEC[iCoord];
		}
		Vector3D normalVector{ globalComponents };
		c_normalVector = std::shared_ptr<Vector3D>(&normalVector);
	}
}
