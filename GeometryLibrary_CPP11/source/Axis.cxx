// baris.albayrak.ieee@gmail.com

#include "Macros.h"
#include "GeometryObject.hxx"
#include "GeometryMath.hxx"
#include "GeometryParameters.hxx"
#include "GeometryException.hxx"
#include "ReferenceObject.hxx"
#include "CoordSystem.hxx"
#include "GlobalCoordSystem.hxx"
#include "Point3D.hxx"
#include "Point2D.hxx"
#include "Point3D.hxx"
#include "Vector3D.hxx"
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
		std::shared_ptr<Point3D> thePassingPoint,
		std::shared_ptr<Vector3D> theDirectionVector)
		:
		GeometryObject(),
		c_passingPoint{ thePassingPoint },
		c_directionVector{ theDirectionVector } { }

	/// <summary>
	/// Ctor
	/// </summary>
	Axis::Axis(
		std::shared_ptr<Point3D> thePoint0,
		ARGCOPY(Point3D) thePoint1)
		: GeometryObject(),
		c_passingPoint{ thePoint0 }
	{
		c_directionVector = std::make_shared<Vector3D>(*thePoint0 , thePoint1);
	}

	/// <summary>
	/// Ctor using the Equation Coefficients
	/// Follows RAII idiom
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Axis::Axis(const std::array<std::array<double, 3>, 2>& theEC)
		: GeometryObject()
	{
		applyEC(theEC);
	}

	/// <summary>
	/// Ctor using the Equation Coefficients
	/// Follows RAII idiom
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Axis::Axis(const std::vector<std::vector<double, std::allocator<double>>>& theEC)
		: GeometryObject()
	{
		applyEC(theEC);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool Axis::operator==(const Axis& rhs) const
	{
		if (&rhs == this) return true;

		auto EC1 = getEC();
		auto EC2 = rhs.getEC();
		if (
			!std::equal(
				EC1[0].cbegin(),
				EC1[0].cend(),
				EC2[0].cbegin(),
				[](double i, double j) {return GeometryMath::equal_g(i, j); })) return false;
		return std::equal(
			EC1[1].cbegin(),
			EC1[1].cend(),
			EC2[1].cbegin(),
			[](double i, double j) {return GeometryMath::equal_g(i, j); });
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Axis::operator!=(const Axis& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincidence.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Axis::operator+=(const Axis& rhs) const
	{
		if (&rhs == this) return true;

		auto EC1 = getEC();
		auto EC2 = rhs.getEC();
		if (
			!std::equal(
				EC1[0].cbegin(),
				EC1[0].cend(),
				EC2[0].cbegin(),
				[](double i, double j) {return GeometryMath::equal_g(i, j); })) return false;
		if (
			std::equal(
				EC1[1].cbegin(),
				EC1[1].cend(),
				EC2[1].cbegin(),
				[](double i, double j) {return GeometryMath::equal_g(i, j); })) return true;
		std::array<double, 3> directionVector2{ {} };
		std::transform(
			EC2[1].cbegin(),
			EC2[1].cend(),
			directionVector2.begin(),
			[](double i) { return -i; });
		return std::equal(
			EC1[1].cbegin(),
			EC1[1].cend(),
			directionVector2.cbegin(),
			[](double i, double j) {return GeometryMath::equal_g(i, j); });
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Axis::operator-=(const Axis& rhs) const
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Axis, Line, Circle and Plane types do not have additional 2D and 3D types (e.g. Circle2D)
	/// Hence, these types do not have is2D and is3D methods.
	/// They are assumed 3D by default.
	/// See project main docstring in GeometryObject.hxx for more detailed description.
	/// However, these types can be geometrically 2D.
	/// This method defines the conditions to have a 2D object wrt the global CS.
	/// </summary>
	bool Axis::is2D() const
	{
		auto EC = getEC();
		if (GeometryMath::zero_g(EC[0][2]) && GeometryMath::zero_g(EC[1][2])) {
			return true;
		}
		return false;
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
	/// See += operator docstring
	/// </summary>
	bool Axis::equalsGeometrically(ARGCOPY(Axis) theAxis) const
	{
		return operator+=(theAxis);
	}

	/// <summary>
	/// Inspect - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::inspectEC(const std::array<std::array<double, 3>, 2>& theEC) const
	{
		if (
			std::all_of(
				theEC[1].cbegin(),
				theEC[1].cend(),
				[](double i) { return GeometryMath::zero_g(i); }))
		{
			throw ZeroVectorException();
		}
	}

	/// <sumrnar y>
	/// Inspect - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::inspectEC(const std::vector<std::vector<double, std::allocator<double>>>& theEC) const
	{
		if (
			std::all_of(
				theEC[1].cbegin(),
				theEC[1].cend(),
				[](double i) { return GeometryMath::zero_g(i); }))
		{
			throw ZeroVectorException();
		}
	}

	/// <summary>
	/// Getter - Passing point
	/// </summary>
	auto Axis::getPassingPoint() const -> std::shared_ptr<Point3D>
	{
		return c_passingPoint;
	}

	/// <summary>
	/// Getter - Direction vector
	/// </summary>
	auto Axis::getDirectionVector() const -> std::shared_ptr<Vector3D>
	{
		return c_directionVector;
	}

	/// <summary>
	/// Getter - Equation Coefficients (EC)
	/// c_EC[0][i] = Passing point coord in ith direction
	/// c_EC[1][i] = Direction vector component in ith direction
	/// </summary>
	auto Axis::getEC() const -> std::array<std::array<double, 3>, 2>
	{
		std::array<std::array<double, 3>, 2> EC = {
			std::array<double, 3>{{}},
			std::array<double, 3>{{}} };
		auto passingPointCoords = c_passingPoint->getGlobalCoords();
		auto directionVectorComponents = c_directionVector->getGlobalComponents();
		std::copy(passingPointCoords.cbegin(), passingPointCoords.cend(), EC[0].begin());
		std::copy(directionVectorComponents.cbegin(), directionVectorComponents.cend(), EC[1].begin());
		return EC;
	}

	/// <summary>
	/// Creates a point by translating the passing point of the axis by the direction vector
	/// </summary>
	auto Axis::createPoint(const double& theFactor) const -> std::shared_ptr<Point3D>
	{
		if (GeometryMath::zero_g(theFactor))
		{
			return c_passingPoint;
		}
		return getDirectionVector()->transformPoint(*getPassingPoint(), theFactor);
	}

	/// <summary>
	/// Returns the location vector of the passing point of the axis.
	/// Location vector is the vector from the global CS origin to the point
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	auto Axis::getPassingPointAsVector() const
	{
		return std::make_shared<Vector3D>(getPassingPoint()->getGlobalCoords());
	}

	/// <summary>
	/// Returns the location vector of the point on the axis defined as:
	///		Transform the passing point using createPoint method.
	///		Location vector is the vector from the global CS origin to the point
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	auto Axis::getPointAsVector(const double& theFactor) const
	{
		return std::make_shared<Vector3D>(createPoint(theFactor)->getGlobalCoords());
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::setMembers(const std::array<std::array<double, 3>, 2>& theEC)
	{
		applyEC(theEC);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Axis::setMembers(const std::vector<std::vector<double, std::allocator<double>>>& theEC)
	{
		applyEC(theEC);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setMembers(
		const std::shared_ptr<Point3D>& thePassingPoint,
		const std::shared_ptr<Vector3D>& theDirectionVector)
	{
		c_passingPoint = thePassingPoint;
		c_directionVector = theDirectionVector;
	}

	/// <summary>
	/// Setter - Passing point
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setPassingPoint(const std::shared_ptr<Point3D>& thePassingPoint)
	{
		c_passingPoint = thePassingPoint;
	}

	/// <summary>
	/// Setter - Direction vector
	/// </summary>
	/// <exception> NullptrException </exception>
	void Axis::setDirectionVector(const std::shared_ptr<Vector3D>& theDirectionVector)
	{
		c_directionVector = theDirectionVector;
	}

	/// <summary>
	/// Shortcut to the corresponding method of Vector3D
	/// </summary>
	bool Axis::isParallel(ARGCOPY(Axis) theAxis) const
	{
		return getDirectionVector()->isParallel(*theAxis.getDirectionVector());
	}

	/// <summary>
	/// Shortcut to the corresponding melhod of Vector3D
	/// </summary>
	bool Axis::isInTheSameDirection(ARGCOPY(Axis) theAxis) const
	{
		return getDirectionVector()->isInTheSameDirection(*theAxis.getDirectionVector());
	}

	/// <summary>
	/// Shortcut to the corresponding method of Vector3D
	/// </summary>
	bool Axis::isNormal(ARGCOPY(Axis) theAxis) const
	{
		return getDirectionVector()->isNormal(*theAxis.getDirectionVector());
	}

	/// <summary>
	/// Returns if the axes are skew
	/// </summary>
	bool Axis::isSkew(ARGCOPY(Axis) theAxis) const
	{
		return intersect(theAxis).first == -1 ? true : false;
	}

	/// <summary>
	/// Returns if the axes are skew
	/// </summary>
	bool Axis::isSkew(ARGCOPY(Line) theLine) const
	{
		auto intersectionResults{ intersect(*theLine.getAxis()) };
		auto intersectionPoint{ intersectionResults.second };
		if (intersectionResults.first == -1) {
			return true;
		}
		if (intersectionResults.first == 0)
		{
			return !theLine.includes(*intersectionPoint) ? true : false;
		}
		return false;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Axis::includes(ARGCOPY(Point3D) thePoint) const
	{
		auto EC = getEC();

		// Calculate the axis parameter at the point 1ocation
		auto globalCoords{ thePoint.getGlobalCoords() };
		double axisParameter{};
		for (int iCoord = 0; iCoord < GeometryParameters::DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::zero_g(EC[1][iCoord])) {
				axisParameter = (globalCoords[iCoord] - EC[0][iCoord]) / EC[1][iCoord];
				break;
			}
		}

		// Calculate coordinate and check if they are equal to the point coordinates
		for (int iCoord = 0; iCoord < GeometryParameters::DIMENSIONS::D3; iCoord++) {
			double globalCoord{ axisParameter * EC[1][iCoord] + EC[0][iCoord] };
			if (!GeometryMath::equal_g(globalCoord, globalCoords[iCoord])) return false;
		}
		return true;
	}

	/// <summary>
	/// Returns if the axis intersects with the input axis
	/// </summary>
	bool Axis::intersects(ARGCOPY(Axis) theAxis) const
	{
		return intersect(theAxis).first == 0;
	}

	/// <summary>
	/// Returns if the axis intersects with the input line
	/// </summary>
	bool Axis::intersects(ARGCOPY(Line) theLine) const
	{
		auto axis1{ theLine.getAxis() };
		auto intersectionResults{ intersect(*axis1) };
		if (intersectionResults.first != 0) return false;
		return theLine.includes(*intersectionResults.second);
	}

	/// <summary>
	/// Returns if the axis coincides with the input axis
	/// </summary>
	bool Axis::coincides(ARGCOPY(Axis) theAxis) const
	{
		return intersect(theAxis).first == 1;
	}

	/// <summary>
	/// Returns if the axis coincides with the input line
	/// </summary>
	bool Axis::coincides(ARGCOPY(Line) theLine) const
	{
		return intersect(*theLine.getAxis()).first == 1;
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input axis
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	auto Axis::intersect(ARGCOPY(Axis) theAxis) const -> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>
	{
		try {
			return intersectBase(
				theAxis,
				getDirectionVector()->crossProduct(*theAxis.getDirectionVector())->getGlobalComponents());
		} catch (ZeroVectorException&) {
			if (includes(*theAxis.getPassingPoint()))
			{
				return std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>{ GeometryParameters::INTERSECTION1::Coincides1, nullptr };
			}
			return std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>{ GeometryParameters::INTERSECTION1::Skew1, nullptr };
		}
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	auto Axis::intersect(ARGCOPY(Line) theLine) const -> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>
	{
		// The axis
		auto intersectionResults{ intersect(*theLine.getAxis()) };
		if (intersectionResults.first != 0)
			return std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>{ intersectionResults.first, nullptr };

		// The line
		intersectionResults.first = (
			theLine.includes(*intersectionResults.second)
			? GeometryParameters::INTERSECTION1::Intersects1
			: GeometryParameters::INTERSECTION1::Skew1);
		return intersectionResults;
	}

	/// <summary>
	/// Projects the input point onto the axis
	/// </summary>
	auto Axis::project(ARGCOPY(Point3D) thePoint) const -> std::shared_ptr<Point3D>
	{
		// Determine the projection
		Point3D passingPoint = *getPassingPoint();
		Vector3D vector0 { passingPoint, thePoint }; // Vector from the passing point to the point
		double distance { vector0.dotProduct(*getDirectionVector()) }; // Distance from the passing point to the projection
		Vector3D vector1 { passingPoint }; // The location vector of the passing point
		std::shared_ptr<Vector3D> vector2 { vector1.add(*getDirectionVector()->multiply(distance)) }; // The location vector of the projection
		return std::make_shared<Point3D>(
				vector2->getReferenceCoordSystem(),
				vector2->getLocalComponents());
	}

	/// <summary>
	/// Calculates distance to an axis
	/// </summary>
	double Axis::calculateDistance(ARGCOPY(Axis) theAxis) const
	{
		// The vector between the passing points
		Vector3D vectorPassingPoint0{ *getPassingPoint() };
		Vector3D vectorPassingPoint1{ *theAxis.getPassingPoint() };
		auto vectorBetweenLines = vectorPassingPoint0.subtruct(vectorPassingPoint1);
		try {
			// The cross product of the direction vectors
			auto crossProduct0 = getDirectionVector()->crossProduct(*theAxis.getDirectionVector());

			// Not parallel: d = (V1 x V2) . (R2 - R1) / |V1 x V2|    V, R1 and R2 are vectors and || demonstrates the magnitude
			// where V is the direction vector, R1 and R2 are the position vectors of the passing points of the lines
			return std::fabs(crossProduct0->dotProduct(*vectorBetweenLines) / crossProduct0->getMagnitude());
		} catch (ZeroVectorException&) {
			return (
				getDirectionVector()->crossProduct(*vectorBetweenLines)->getMagnitude() /
				getDirectionVector()->getMagnitude());
		}
	}

	/// <summary>
	/// Calculates distance to a point.
	/// </summary>
	double Axis::calculateDistance(ARGCOPY(Point3D) thePoint) const
	{
		return thePoint.calculateDistance(*project(thePoint));
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
		auto axis1{ theLine.getAxis() };
		auto distanceInfinite{ calculateDistance(*axis1) };

		// The axes are parallel
		if (isParallel(*axis1)) return distanceInfinite;

		// The axes intersect
		if (GeometryMath::zero_g(distanceInfinite)) {
			auto intersectionResults{ intersect(*axis1) };
			if (theLine.includes(*intersectionResults.second)) return 0.;

			return std::fmin(
				calculateDistance(*theLine.getEndPoint0()),
				calculateDistance(*theLine.getEndPoint1()));
		}

		// The axes are skew
		auto closestPoints{ findClosestPoints(*axis1) };
		if (theLine.includes(*closestPoints[1]))
		{
			return closestPoints[0]->calculateDistance(*closestPoints[1]);
		}
		return std::fmin(
			closestPoints[0]->calculateDistance(*theLine.getEndPoint0()),
			closestPoints[0]->calculateDistance(*theLine.getEndPoint1()));
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
	auto Axis::findClosestPoints(ARGCOPY(Axis) theAxis) const -> std::vector<std::shared_ptr<Point3D>>
	{
		// Initialize the output
		std::vector<std::shared_ptr<Point3D>> outPoints;

		// Get the passing points and the direction vectors
		auto passingPoint0 = getPassingPoint();
		auto passingPoint1 = theAxis.getPassingPoint();
		auto directionVector0{ getDirectionVector() };
		auto directionVector1{ theAxis.getDirectionVector() };
		try {
			auto crossProduct0 = directionVector1->crossProduct(*directionVector0);
		} catch (ZeroVectorException&) {
			// Two lines are parallel:
			// Create a line perpandicular to both lines at the passing point of this
			// The passing point and the intersection of the new line with the input line are the outputs
			outPoints.push_back(getPassingPoint());
			auto dummy = Vector3D(*passingPoint0, *passingPoint1);
			auto crossProduct1 = dummy.crossProduct(*directionVector0);
			auto crossProduct2 = directionVector0->crossProduct(*crossProduct1);
			Axis axis { getPassingPoint(), crossProduct2 };
			auto intersectionResults{ intersect(axis)};
			outPoints.push_back(intersectionResults.second);

			return outPoints;
		}

		// Skew lines: P0 + t0V0 + t2V2 = P1 + tM (see the method docstring)
		// Get the coords of the passing points and the components of the direction vectors
		auto crossProduct0 = directionVector1->crossProduct(*directionVector0);
		auto P0{ passingPoint0->getGlobalCoords() };
		auto P1{ passingPoint1->getGlobalCoords() };
		auto V0{ directionVector0->getGlobalComponents() };
		auto V1{ directionVector1->getGlobalComponents() };
		auto V2{ crossProduct0->getGlobalComponents() };

		// Store the coefficients of the thre equations in a matrix for the solution
		std::array<std::array<double, 3>, 3> coefficients = {
			std::array<double, 3>{{}},
			std::array<double, 3>{{}},
			std::array<double, 3>{{}} };
		std::array<double, 3> resultants = { {} };
		for (int iCoord = 0; iCoord < GeometryParameters::DIMENSIONS::D3; iCoord++) {
			coefficients[iCoord] = { V0[iCoord], V1[iCoord] * (-1), V2[iCoord] };
			resultants[iCoord] = P1[iCoord] - P0[iCoord];
		}

		// Solve the equations
		auto inverseCoefficients{ GeometryMath::calculateMatrixInverseS33(coefficients) };
		auto roots{
			GeometryMath::multiplyMatrixToVector<double, 3>(
				inverseCoefficients,
				resultants) };
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
		auto EC = getEC();
		if (GeometryMath::zero_g(EC[1][theIndexCoord1])) throw AssymptoticLineException();
		return (theCoord - EC[0][theIndexCoord1]) * EC[1][theIndexCoord0] / EC[1][theIndexCoord1] + EC[0][theIndexCoord0];
	}

	/// <summary>
	/// Calculate the coord of the point on the line - Input Y, Requested X
	/// Throws exception if the line is parallel to the requested axis or constant in the input axis
	/// </summary>
	/// <exception> AssymptoticLineException </exception>
	double Axis::calculateCoordX_fromCoordY(const double& theCoordY) const
	{
		if (getDirectionVector()->isParallel(*Vector3D::createUnitVectorX()))
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
		if (getDirectionVector()->is2D())
		{
			throw DimensionalityException();
		}

		if (getDirectionVector()->isParallel(*Vector3D::createUnitVectorX()))
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
		if (getDirectionVector()->isParallel(*Vector3D::createUnitVectorY()))
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
		if (getDirectionVector()->is2D())
		{
			throw DimensionalityException();
		}

		if (getDirectionVector()->isParallel(*Vector3D::createUnitVectorY()))
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
		if (getDirectionVector()->is2D()) throw DimensionalityException();

		if (getDirectionVector()->isParallel(*Vector3D::createUnitVectorZ()))
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
		if (getDirectionVector()->is2D()) throw DimensionalityException();

		if (getDirectionVector()->isParallel(*Vector3D::createUnitVectorZ()))
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
	auto Axis::intersectBase(
		ARGCOPY(Axis) theAxis,
		const std::array<double, 3>& theCrossProductComponents) const
		-> std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>
	{
		auto EC1 = getEC();

		// Determine the axis parameters at the intersection
		// No need to inspect exception as the parallelism is inspected before
		std::array<std::array<double, 3>, 2> EC2{ theAxis.getEC() };
		double p0_x{ EC1[0][0] };
		double p0_y{ EC1[0][1] };
		double p0_z{ EC1[0][2] };
		double p1_x{ EC2[0][0] };
		double p1_y{ EC2[0][1] };
		double p1_z{ EC2[0][2] };
		double V0_x{ EC1[1][0] };
		double V0_y{ EC1[1][1] };
		double V0_z{ EC1[1][2] };
		double V1_x{ EC2[1][0] };
		double V1_y{ EC2[1][1] };
		double V1_z{ EC2[1][2] };
		double axisParameter0;
		double axisParameter1;
		if (!GeometryMath::zero_g(theCrossProductComponents[0])) {
			axisParameter0 = (V1_z * (p1_y - p0_y) - V1_y * (p1_z - p0_z)) / (V0_y * V1_z - V0_z * V1_y);
			axisParameter1 = (V0_z * (p1_y - p0_y) - V0_y * (p1_z - p0_z)) / (V0_y * V1_z - V0_z * V1_y);
		}
		else if (!GeometryMath::zero_g(theCrossProductComponents[1])) {
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
			!GeometryMath::equal_g(globalCoord0_x, globalCoord1_x) ||
			!GeometryMath::equal_g(globalCoord0_y, globalCoord1_y) ||
			!GeometryMath::equal_g(globalCoord0_z, globalCoord1_z))
		{
			return std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>>{ GeometryParameters::INTERSECTION1::Skew1, nullptr };
		}

		std::array<double, 3> intersectionCoords {{ globalCoord0_x, globalCoord0_y, globalCoord0_z}};
		auto intersection = std::make_shared<Point3D>(intersectionCoords);

		// Outputs
		auto referenceCoordSystem0 = getDirectionVector()->getReferenceCoordSystem();
		auto referenceCoordSystem1 = theAxis.getDirectionVector()->getReferenceCoordSystem();
		std::pair<GeometryParameters::INTERSECTION1, std::shared_ptr<Point3D>> outIntersectionResults;
		outIntersectionResults.first = GeometryParameters::INTERSECTION1::Intersects1;
		if (is2D() != theAxis.is2D() || referenceCoordSystem0 != referenceCoordSystem1)
		{
			outIntersectionResults.second = intersection;
		}
		else
		{
			outIntersectionResults.second = std::make_shared<Point3D>(
				referenceCoordSystem0,
				referenceCoordSystem0->measurePointCoords(*intersection));
		}
		return outIntersectionResults;
	}

	/// <summary>
	/// Determine the passing point and direction vector corresponding to the input EC
	/// Set the members accordingly
	/// </summary>
	void Axis::applyEC(const std::array<std::array<double, 3>, 2>& theEC)
	{
		inspectEC(theEC);

		std::array<double, 3> passingPointCoords = { {} };
		std::array<double, 3> directionVectorComponents = { {} };
		copy(theEC[0].cbegin(), theEC[0].cend(), passingPointCoords.begin());
		copy(theEC[1].cbegin(), theEC[1].cend(), directionVectorComponents.begin());

		c_passingPoint = std::make_shared<Point3D>(passingPointCoords);
		c_directionVector = std::make_shared<Vector3D>(directionVectorComponents);
	}

	void Axis::applyEC(const std::vector<std::vector<double, std::allocator<double>>>& theEC)
	{
		inspectEC(theEC);

		std::array<double, 3> passingPointCoords = { {} };
		std::array<double, 3> directionVectorComponents = { {} };
		copy(theEC[0].cbegin(), theEC[0].cend(), passingPointCoords.begin());
		copy(theEC[1].cbegin(), theEC[1].cend(), directionVectorComponents.begin());

		c_passingPoint = nullptr;
		c_passingPoint = std::make_shared<Point3D>(passingPointCoords);
		c_directionVector = nullptr;
		c_directionVector = std::make_shared<Vector3D>(directionVectorComponents);
	}
}
