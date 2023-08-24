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
	/// Ctor using the Equation Coefficients
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Plane::Plane(const arrayS4& theEC)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(theEC);
	}

	/// <summary>
	/// Ctor using the Equation Coefficients
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Plane::Plane(const vectorInput1D& theEC)
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
		PointBase& thePoint,
		VectorBase& theNormalVector)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(thePoint, theNormalVector);
	}

	/// <summary>
	/// Ctor using the passing point and two in-plane vectors
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	Plane::Plane(
		PointBase& thePoint,
		ARGCOPY(VectorBase) theInPlaneVector0,
		ARGCOPY(VectorBase) theInPlaneVector1)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		Vector3D vector { theInPlaneVector0.crossProduct(theInPlaneVector1) };
		c_passingPoint = shared_ptr<PointBase>(&thePoint);
		c_normalVector = shared_ptr<VectorBase>(&vector);
		//setMembers(thePoint, theInPlaneVector0.crossProduct(theInPlaneVector1));
	}

	/// <summary>
	/// Ctor using three points
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> CoincidenceException </exception>
	/// <exception> ColinearPointsException </exception>
	Plane::Plane(
		PointBase& thePoint0,
		ARGCOPY(PointBase) thePoint1,
		ARGCOPY(PointBase) thePoint2)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		// Inspect inputs
		if (thePoint0.coincides(thePoint1)) throw CoincidenceException();
		if (thePoint0.coincides(thePoint2)) throw CoincidenceException();

		Axis axis { thePoint0, thePoint1 };
		if (axis.includes(thePoint2)) {
			throw ColinearPointsException();
		}

		// Create two vectors in the plane
		Vector3D vector0{ thePoint0, thePoint1 };
		Vector3D vector1{ thePoint0, thePoint2 };

		Vector3D vector { vector0.crossProduct(vector1) };
		c_passingPoint = shared_ptr<PointBase>(&thePoint0);
		c_normalVector = shared_ptr<VectorBase>(&vector);
		//setMembers(thePoint0, vector0.crossProduct(vector1));
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Plane::Plane(const Plane& rhs)
	{
		copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Plane& Plane::operator=(const Plane& rhs) {
		if (&rhs == this) return *this;

		copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Plane::Plane(Plane&& rhs) noexcept
	{
		copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Plane& Plane::operator=(Plane&& rhs) noexcept
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
	bool Plane::operator==(const Plane& rhs)
	{
		if (&rhs == this) return true;
		if (*c_passingPoint != rhs.getPassingPoint()) return false;
		return *c_normalVector == rhs.getNormalVector();
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Plane::operator!=(const Plane& rhs)
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Plane::operator+=(const Plane& rhs)
	{
		if (&rhs == this) return true;
		if (!c_normalVector->equalsGeometrically(rhs.getNormalVector())) return false;
		return includes(rhs.getPassingPoint());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Plane::operator-=(const Plane& rhs)
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	Plane::~Plane()
	{
		Destroy();
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void Plane::copyBase(const Plane& rhs)
	{
		GeometryObject::copyBase(rhs);
		std::copy(std::begin(rhs.getEC()), std::end(rhs.getEC()), std::begin(c_EC));

		PointBase passingPoint = rhs.getPassingPoint();
		VectorBase normalVector { 
			rhs.getNormalVector().getDimensionCount(),
			rhs.getNormalVector().getReferenceCoordSystem(),
			rhs.getNormalVector().getLocalComponents() };
		c_passingPoint = shared_ptr<PointBase>(&passingPoint);
		c_normalVector = shared_ptr<VectorBase>(&normalVector);
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
			!GeometryMath::equals(c_passingPoint->getLocalCoordZ(), 0., c_toleranceGeneral))
		{
			return false;
		}
		if (c_passingPoint->getReferenceCoordSystem().getAxisAsVectorZ().isParallel(*c_normalVector))
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
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	void Plane::Destroy()
	{
		c_passingPoint = nullptr;
		c_normalVector = nullptr;
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool Plane::equalsBase(ARGCOPY(Plane) thePlane) const
	{
		return true;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool Plane::equals(ARGCOPY(Plane) thePlane) const
	{
		if (this == &thePlane) return true;
		if (!equalsBase(thePlane)) return false;
		if (!c_passingPoint->equals(thePlane.getPassingPoint())) return false;
		if (!c_normalVector->equals(thePlane.getNormalVector())) return false;
		return true;
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool Plane::equalsGeometrically(ARGCOPY(Plane) thePlane) const
	{
		if (this == &thePlane) return true;
		if (!equalsBase(thePlane)) return false;

		// The passing points may be different even for eqivalent planes.
		// Hence, unequivalence of passing points does not mean unequivalence of the planes.
		// However, equivalence of passing points and normal vectors means equivalent planes.
		if (c_normalVector->equalsGeometrically(thePlane.getNormalVector())) return false;
		return includes(thePlane.getPassingPoint());
	}

	/// <summary>
	/// Inspect Equation Coefficients (EC)
	/// At least one of the 1st three coefficients shall be non-zero.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::inspectEC(const arrayS4& theEC)
	{
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::equals(theEC[iCoord], 0., c_toleranceGeneral)) {
				return;
			}
		}
		throw ZeroVectorException();
	}

	/// <summary>
	/// Inspect Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::inspectEC(const vectorInput1D& theEC)
	{
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			if (!GeometryMath::equals(theEC[iCoord], 0., c_toleranceGeneral)) {
				return;
			}
		}
		throw ZeroVectorException();
	}

	/// </summary>
	/// Getter - Passing point
	/// </summary>
	PointBase& Plane::getPassingPoint() const
	{
		return *c_passingPoint;
	}

	/// <summary>
	/// Getter - Normal vector
	/// </summary>
	VectorBase& Plane::getNormalVector() const
	{
		return *c_normalVector;
	}

	/// <summary>
	/// Getter - EC
	/// </summary>
	arrayS4 Plane::getEC() const
	{
		return c_EC;
	}

	/// <summary>
	/// Returns the reference CS which is common to all ReferenceObject members.
	///		Returns a handle with nullptr if the ReferenceObject members have different reference CSs.
	/// See module docstring in Axis.hxx for the details.
	/// </summary>
	CoordSystem* Plane::getCommonReferenceCoordSystem() const
	{
		if (c_passingPoint->getReferenceCoordSystem() == c_normalVector->getReferenceCoordSystem())
		{
			return &(c_passingPoint->getReferenceCoordSystem());
		}
		return nullptr;
	}

	/// <summary>
	/// Setter for the members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::setMembers(const arrayS4& theEC)
	{
		setEC(theEC);
	}

	/// <summary>
	/// Setter for the members
	/// Protected method used in this class and child classes only.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::setMembers(const vectorInput1D& theEC)
	{
		setEC(theEC);
	}

	/// <summary>
	/// Setter for the members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void Plane::setMembers(
		PointBase& thePassingPoint,
		VectorBase& theNormalVector)
	{
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		c_normalVector = std::shared_ptr<VectorBase>(&theNormalVector);
		updateEC();
	}

	/// <summary>
	/// Setter - Passing point
	/// Creates a point at the globalorigin if the input point is null.
	/// </summary>
	void Plane::setPassingPoint(PointBase& thePassingPoint)
	{
		c_passingPoint = std::shared_ptr<PointBase>(&thePassingPoint);
		updateEC();
	}

	/// <summary>
	/// Setter - Direction vector
	/// Creates a unit vector (z = 1} at the globalorigin if the input vector is null.
	/// </summary>
	void Plane::setNormalVector(VectorBase& theNormalVector)
	{
		c_normalVector = std::shared_ptr<VectorBase>(&theNormalVector);
		updateEC();
	}

	/// <sumrnary>
	/// Setter - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::setEC(const arrayS4& theEC)
	{
		inspectEC(theEC);
		c_EC = theEC;
		applyEC();
	}

	/// <summary>
	/// Setter - Equation Coefficients (EC)
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	void Plane::setEC(const vectorInput1D& theEC)
	{
		inspectEC(theEC);
		c_EC = GeometryMath::convertVectorToArray1DS4(theEC);
		applyEC();
	}

	/// <summary>
	/// Returns if the input intersets
	/// </summary>
	bool Plane::intersects(ARGCOPY(Axis) theAxis) const
	{
		std::pair<int, Point3D*> intersectionResults { intersect(theAxis) };
		return bool(intersectionResults.first == 0);
	}

	/// <summary>
	/// Returns if the input intersets
	/// </summary>
	bool Plane::intersects(ARGCOPY(Line) theLine) const
	{
		std::pair<int, Point3D*> intersectionResults { intersect(theLine) };
		bool outResult{ intersectionResults.first == 0 ? true : false };
		return outResult;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Plane::includes(ARGCOPY(PointBase) thePoint) const
	{
		// Calculate z-coordinates of the point on the plane having x and y coordinates
		// same as the input point and compare with the point's z-coordinate
		arrayS3 globalCoords{ thePoint.getGlobalCoords() };
		double coordX{ globalCoords[0] };
		double coordY{ globalCoords[1] };
		std::pair<bool, double> coordPair{ calculateCoordZ(coordX, coordY) };
		if (
			!coordPair.first ||
			GeometryMath::equals(
				coordPair.second,
				globalCoords[2],
				c_toleranceGeneral))
		{
			return true;
		}
		return false;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Plane::includes(ARGCOPY(Axis) theAxis) const
	{
		// Check if passing point of the line is included
		if (!includes(theAxis.getPassingPoint())) return false;

		// Check if the plane normal vector and the line direction vector are normal to each other
		if (c_normalVector->isNormal(theAxis.getDirectionVector())) return true;
		return false;
	}

	/// <summary>
	/// Returns if the input is included
	/// </summary>
	bool Plane::includes(ARGCOPY(Line) theLine) const
	{
		// Check if passing point of the line is included
		if (!includes(theLine.getEndPoint0())) return false;

		// Check if the plane normal vector and the line direction vector are normal to each other
		if (c_normalVector->isNormal(theLine.getDirectionVector())) return true;
		return false;
	}

	/// <summary>
	/// Calculates the coord X
	/// </summary>
	std::pair<bool, double> Plane::calculateCoordX(const double& theCoordY, const double& theCoordZ) const
	{
		// Initialize outputs
		std::pair<bool, double> outResults{ false, 0. };

		// Parameterize inputs
		double coord1{ theCoordY };
		double coord2{ theCoordZ };
		double EC0{ c_EC[0] };
		double EC1{ c_EC[1] };
		double EC2{ c_EC[2] };
		double EC_3{ c_EC[3] };

		// Check if the plane includes the input coordinates
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral)) return outResults;
		outResults.first = true;

		// The output coordinate
		outResults.second = -(EC1 * coord1 + EC2 * coord2 + EC_3) / EC0;
		return outResults;
	}

	/// <summary>
	/// Calculates the coord Y
	/// </summary>
	std::pair<bool, double> Plane::calculateCoordY(const double& theCoordZ, const double& theCoordX) const
	{
		// Initialize outputs
		std::pair<bool, double> outResults{ false, 0. };

		// Parameterize inputs
		double coord1{ theCoordZ };
		double coord2{ theCoordX };
		double EC0{ c_EC[1] };
		double EC1{ c_EC[2] };
		double EC2{ c_EC[0] };
		double EC_3{ c_EC[3] };

		// Check if the plane includes the input coordinates
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral)) return outResults;
		outResults.first = true;

		// The output coordinate
		outResults.second = -(EC1 * coord1 + EC2 * coord2 + EC_3) / EC0;
		return outResults;
	}

	/// <summary>
	/// Calculates the coord Z
	/// </summary>
	std::pair<bool, double> Plane::calculateCoordZ(const double& theCoordX, const double& theCoordY) const
	{
		// Initialize outputs
		std::pair<bool, double> outResults{ false, 0. };

		// Parameterize inputs
		double coord1{ theCoordX };
		double coord2{ theCoordY };
		double EC0{ c_EC[2] };
		double EC1{ c_EC[0] };
		double EC2{ c_EC[1] };
		double EC_3{ c_EC[3] };

		// Check if the plane includes the input coordinates
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral)) return outResults;
		outResults.first = true;

		// The output coordinate
		outResults.second = -(EC1 * coord1 + EC2 * coord2 + EC_3) / EC0;
		return outResults;
	}

	/// <summary>
	/// Calculates the coords X and Y
	/// </summary>
	std::pair<bool, std::array<double, 2>> Plane::calculateCoordXY(const double& theCoordZ) const
	{
		// Initialize outputs
		std::pair<bool, std::array<double, 2>> outResults{ false, {0., 0.} };

		// Parameterize inputs
		double coord2{ theCoordZ };
		double EC0{ c_EC[0] };
		double EC1{ c_EC[1] };
		double EC2{ c_EC[2] };
		double EC_3{ c_EC[3] };

		// Check if the plane includes the input coordinate
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral) && GeometryMath::equals(EC1, 0., c_toleranceGeneral))
		{
			return outResults;
		}
		outResults.first = true;

		// Check if the plane is dependent only on the st output direction
		if (GeometryMath::equals(EC1, 0., c_toleranceGeneral) && GeometryMath::equals(EC2, 0., c_toleranceGeneral))
		{
			outResults.second[0] = -EC_3 / EC0;
			return outResults;
		}

		// Check if the plane is dependent only on the 2nd output direction
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral) && GeometryMath::equals(EC2, 0., c_toleranceGeneral))
		{
			outResults.second[1] = -EC_3 / EC1;
			return outResults;
		}

		// Check if the plane is independent of the st output direction
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral)) {
			outResults.second[1] = -(EC_3 + EC2 * coord2) / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 2nd output direction
		if (GeometryMath::equals(EC1, 0., c_toleranceGeneral)) {
			outResults.second[0] = -(EC_3 + EC2 * coord2) / EC0;
			return outResults;
		}

		// All exceptional cases inspected.
		// All directions are dependent
		// Keep the 1st output direction with zero initial value and calculate the 2nd direction correspondingly
		outResults.second[1] = -(EC_3 + EC2 * coord2) / EC1;
		return outResults;
	}

	/// <summary>
	/// Calculates the coords Y and Z
	/// </summary>
	std::pair<bool, std::array<double, 2>> Plane::calculateCoordYZ(const double& theCoordX) const
	{
		// Initialize outputs
		std::pair<bool, std::array<double, 2>> outResults{ false, {0., 0.} };

		// Parameteri ze inputs
		double coord2{ theCoordX };
		double EC0{ c_EC[1] };
		double EC1{ c_EC[2] };
		double EC2{ c_EC[0] };
		double EC_3{ c_EC[3] };

		// Check if the plane includes the input coordinate
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral) && GeometryMath::equals(EC1, 0., c_toleranceGeneral))
		{
			return outResults;
		}
		outResults.first = true;

		// Check if the plane is dependent only on the st output direction
		if (GeometryMath::equals(EC1, 0., c_toleranceGeneral) && GeometryMath::equals(EC2, 0., c_toleranceGeneral))
		{
			outResults.second[0] = -EC_3 / EC0;
			return outResults;
		}

		// Check if the plane is dependent only on the 2nd output direction
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral) && GeometryMath::equals(EC2, 0., c_toleranceGeneral))
		{
			outResults.second[1] = -EC_3 / EC1;
			return outResults;
		}

		// Check if the plane is independent of the st output direction
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral)) {
			outResults.second[1] = -(EC_3 + EC2 * coord2) / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 2nd output direction
		if (GeometryMath::equals(EC1, 0., c_toleranceGeneral)) {
			outResults.second[0] = -(EC_3 + EC2 * coord2) / EC0;
			return outResults;
		}

		// All exceptional cases inspected.
		// All directions are dependent
		// Keep the st output direction with zero initial value and calculatethe 2nd direction correspondingly
		outResults.second[1] = -(EC_3 + EC2 * coord2) / EC1;
		return outResults;
	}

	/// <summary>
	/// Calculates the coords Z and X
	/// </summary>
	std::pair<bool, std::array<double, 2>> Plane::calculateCoordZX(const double& theCoordY) const
	{
		// Initialize outputs
		std::pair<bool, std::array<double, 2>> outResults{ false, {0., 0.} };

		// Parameterize inputs
		double coord2{ theCoordY };
		double EC0{ c_EC[2] };
		double EC1{ c_EC[0] };
		double EC2{ c_EC[1] };
		double EC_3{ c_EC[3] };

		// Check if the plane includes the input coordinate
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral) && GeometryMath::equals(EC1, 0., c_toleranceGeneral))
		{
			return outResults;
		}
		outResults.first = true;

		// Check if the plane is dependent only on the st output direction
		if (GeometryMath::equals(EC1, 0., c_toleranceGeneral) && GeometryMath::equals(EC2, 0., c_toleranceGeneral))
		{
			outResults.second[0] = -EC_3 / EC0;
			return outResults;
		}

		// Check if the plane is dependent only on the 2nd output direction
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral) && GeometryMath::equals(EC2, 0., c_toleranceGeneral))
		{
			outResults.second[1] = -EC_3 / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 1st output direction
		if (GeometryMath::equals(EC0, 0., c_toleranceGeneral)) {
			outResults.second[1] = -(EC_3 + EC2 * coord2) / EC1;
			return outResults;
		}

		// Check if the plane is independent of the 2nd output direction
		if (GeometryMath::equals(EC1, 0., c_toleranceGeneral)) {
			outResults.second[0] = -(EC_3 + EC2 * coord2) / EC0;
			return outResults;
		}

		// All exceptional cases inspected.
		// All directions are dependent
		// Keep the 1st output direction with zero initial value and calculatethe 2nd direction correspondingly
		outResults.second[1] = -(EC_3 + EC2 * coord2) / EC1;
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
	double Plane::calculateDistance(ARGCOPY(PointBase) thePoint) const
	{
		Point3D point { thePoint.getReferenceCoordSystem(), thePoint.getLocalCoords() };
		double value0{};
		double value1{};
		arrayS3 globalCoords{ point.getGlobalCoords() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			value0 += c_EC[iCoord] * globalCoords[iCoord];
			value1 += std::pow(c_EC[iCoord], 2.);
		}

		value0 += c_EC[3];
		value1 = std::pow(value1, 0.5);
		return std::fabs(value0 / value1);
	}

	/// <summary>
	/// Calculates distance to the input axis
	/// </summary>
	double Plane::calculateDistance(ARGCOPY(Axis) theAxis) const
	{
		// Check if the two intersects
		if (!c_normalVector->isNormal(theAxis.getDirectionVector())) return 0.;

		// Find distance from the line's passing point to the plane
		return calculateDistance(theAxis.getPassingPoint());
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
			calculateDistance(theLine.getEndPoint0()), 
			calculateDistance(theLine.getEndPoint1()));
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
	std::pair<INTERSECTION2, Axis*> Plane::intersect(ARGCOPY(Plane) thePlane) const
	{
		// Check parallelism
		if (c_normalVector->isParallel(thePlane.getNormalVector()))
		{
			if (GeometryMath::equals(calculateDistance(*thePlane.c_passingPoint), 0., c_toleranceGeneral))
			{
				return std::pair<INTERSECTION2, Axis*>{ INTERSECTION2::Includes2, nullptr };
			}
			return std::pair<INTERSECTION2, Axis*>{ INTERSECTION2::Skew2, nullptr };
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
		bool exists{};
		double coordX{};
		double coordY{};
		double coordZ{};
		double divisor{};
		double denominator{};

		// Trial 1: x = o.
		coordX = 0.;
		denominator = c_EC[1] * thePlane.c_EC[2] - thePlane.c_EC[1] * c_EC[2];
		if (!GeometryMath::equals(denominator, 0., c_toleranceGeneral)) {
			divisor = c_EC[2] * thePlane.c_EC[3] - thePlane.c_EC[2] * c_EC[3];
			coordY = divisor / denominator;
			std::pair<bool, double> coordPair = thePlane.calculateCoordZ(coordX, coordY);
			exists = coordPair.first;
		}
		else {
			// Trial 2: y = 0.
			coordY = 0.;
			denominator = c_EC[2] * thePlane.c_EC[0] - thePlane.c_EC[2] * c_EC[0];
			if (!GeometryMath::equals(denominator, 0., c_toleranceGeneral)) {
				divisor = c_EC[0] * thePlane.c_EC[3] - thePlane.c_EC[0] * c_EC[3];
				coordZ = divisor / denominator;
				std::pair<bool, double> coordPair = thePlane.calculateCoordX(coordY, coordZ);
				exists = coordPair.first;
			}
			else {
				// Trial 3: z = 0.
				coordZ = 0.;
				denominator = c_EC[0] * thePlane.c_EC[1] - thePlane.c_EC[0] * c_EC[1];
				if (!GeometryMath::equals(denominator, 0., c_toleranceGeneral)) {
					divisor = c_EC[1] * thePlane.c_EC[3] - thePlane.c_EC[1] * c_EC[3];
					coordX = divisor / denominator;
					std::pair<bool, double> coordPair = thePlane.calculateCoordY(coordZ, coordX);
					exists = coordPair.first;
				}
			}
		}
		Point3D passingPoint { arrayS3{ coordX, coordY, coordZ } };

		// The direction vector of the intersection line
		// No need to inspect exception as the parallelism is checl<ed at the beginning of the method
		Vector3D directionVector{ c_normalVector->crossProduct(*thePlane.c_normalVector) };

		// intersection line
		return std::pair<INTERSECTION2, Axis*>{
			INTERSECTION2::Intersects2,
			new Axis(
				passingPoint,
				directionVector) };
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
	std::pair<INTERSECTION2, Point3D*> Plane::intersect(ARGCOPY(Axis) theAxis) const
	{
		// Check parallelism
		double angle{ std::fabs(c_normalVector->calculateAngle(theAxis.getDirectionVector())) };
		if (
			GeometryMath::equals(angle, M_PI / 2., c_toleranceSensitive) ||
			GeometryMath::equals(angle, M_PI * 3. / 2., c_toleranceSensitive))
		{
			if (GeometryMath::equals(calculateDistance(theAxis.getPassingPoint()), 0., c_toleranceGeneral))
			{
				return std::pair<INTERSECTION2, Point3D*>{ INTERSECTION2::Includes2, nullptr };
			}
			return std::pair<INTERSECTION2, Point3D*>{ INTERSECTION2::Skew2, nullptr };
		}

		// Solve the plane's scalar equation (Ax + By + Cz + D = 0) and line's parametric equation
		// (x - x0/Vx = y - y0/Vy = z - z0/Vz = t) together to obtain line's parameter (t)
		// A * (x0 + Vx * t) + 8 * (y0 + Vy * t) + C * (z0 + Vz * t) + D = 0 where t is the only unknown
		// No need to inspect exception as the parallelism is inspected before
		arrayS32 EC1{ theAxis.getEC() };
		double axisParameter{ -(
			(
				c_EC[0] * EC1[0][0] +
				c_EC[1] * EC1[1][0] +
				c_EC[2] * EC1[2][0] +
				c_EC[3]) /
			(
				c_EC[0] * EC1[0][1] +
				c_EC[0] * EC1[1][1] +
				c_EC[0] * EC1[2][1])) };

		// Determine point coords using line1 s parameteric equation with the parameter
		arrayS3 pointCoords{
			EC1[0][0] + EC1[0][1] * axisParameter,
			EC1[1][0] + EC1[1][1] * axisParameter,
			EC1[2][0] + EC1[2][1] * axisParameter };

		// Create point
		return std::pair<INTERSECTION2, Point3D*>{ INTERSECTION2::Intersects2, new Point3D(pointCoords) };
	}

	/// <summary>
	/// Returns the intersection status and the intersection object (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	std::pair<INTERSECTION2, Point3D*> Plane::intersect(ARGCOPY(Line) theLine) const
	{
		// Inspect intersection for the infinite tine corresponding to the line
		std::pair<INTERSECTION2, Point3D*> intersectionResults{ intersect(theLine.getAxis()) };
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
	Point3D Plane::project(ARGCOPY(PointBase) thePoint) const
	{
		// Check if the point is on the plane
		if (includes(thePoint)) {
			return Point3D(thePoint.getReferenceCoordSystem(), thePoint.getLocalCoords());
		}

		// Create a line passing through the input point and in the direction parallel to the planes normal vector
		PointBase point { thePoint };
		Axis axis { point, *c_normalVector };

		// The projection is the intersection point
		std::pair<INTERSECTION2, Point3D*> intersectionResults{ intersect(axis) };
		return *intersectionResults.second;
	}

	/// <summary>
	/// Projects the input vector onto the plane
	/// </summary>
	Vector3D* Plane::project(ARGCOPY(VectorBase) theVector) const
	{
		// Create an arbitray axis using the vector
		Point3D arbitraryPoint { arrayS3{} };
		Point3D transformPoint { theVector.transformPoint(arbitraryPoint, 2.) };
		Axis axis { arbitraryPoint, transformPoint };

		// Project the axis onto the plane
		Axis* projection { project(axis) };
		if (!projection) return nullptr;

		// The direction vector of the projection is the projected vector
		return (Vector3D*)(&projection->getDirectionVector());
	}

	/// <summary>
	/// Projects the input axis onto the plane
	/// </summary>
	Axis* Plane::project(ARGCOPY(Axis) theAxis) const
	{
		// Take two points on the line and project onto the plane
		Point3D pointOnPlane0 { project(theAxis.getPassingPoint()) };
		PointBase pointOnLine0 { theAxis.createPoint(2.) };
		Point3D pointOnLine1 { pointOnLine0.getReferenceCoordSystem(), pointOnLine0.getLocalCoords() };
		Point3D pointOnPlane1 { project(pointOnLine1) };
		if (pointOnPlane0.coincides(pointOnPlane1))
		{
			return nullptr;
		}
		return new Axis(pointOnPlane0, pointOnPlane1);
	}

	/// <summary>
	/// Projects the input line onto the plane
	/// </summary>
	Line* Plane::project(ARGCOPY(Line) theLine) const
	{
		// Take two points on the line and project onto the plane
		Point3D pointOnPlane0 { project(theLine.getEndPoint0()) };
		Point3D pointOnPlane1 { project(theLine.getEndPoint1()) };
		if (pointOnPlane0.coincides(pointOnPlane1))
		{
			return nullptr;
		}

		// Create line using the projection points
		return new Line(pointOnPlane0, pointOnPlane1);
	}

	/// <summary>
	/// Creates a point using Equation Coefficients (EC)
	/// </summary>
	/// <exception> UncaughtException </exception>
	Point3D Plane::createPoint(const double& theFactor) const
	{
		// Set a value to one of the thre coordinates and determine the other two coordinates
		double coordX{};
		bool exists{};
		double coordY{};
		double coordZ{};
		arrayS3 globalCoords{ c_passingPoint->getGlobalCoords() };
		if (!GeometryMath::equals(c_EC[0], 0., c_toleranceGeneral)) {
			coordX = globalCoords[0] + theFactor;
			std::pair<bool, std::array<double, 2>> coordPair{ calculateCoordYZ(coordX) };
			exists = coordPair.first;
		}
		else if (!GeometryMath::equals(c_EC[1], 0., c_toleranceGeneral)) {
			coordY = globalCoords[1] + theFactor;
			std::pair<bool, std::array<double, 2>> coordPair{ calculateCoordZX(coordY) };
			exists = coordPair.first;
		}
		else if (!GeometryMath::equals(c_EC[2], 0., c_toleranceGeneral)) {
			coordZ = globalCoords[2] + theFactor;
			std::pair<bool, std::array<double, 2>> coordPair{ calculateCoordXY(coordZ) };
			exists = coordPair.first;
		}
		else {
			coordX = globalCoords[0] + theFactor;
			std::pair<bool, std::array<double, 2>> coordPair{ calculateCoordYZ(coordX) };
			exists = coordPair.first;
		}
		if (!exists) throw UncaughtException();

		return Point3D (arrayS3{ coordX, coordY, coordZ});
	}

	/// <summary>
	/// Internal method to calculate ECs
	/// Requires non-null value for the passing point and the direction vector members
	/// </summary>
	void Plane::updateEC()
	{
		c_EC[3] = 0.;

		double EC_3 = 0.;
		arrayS3 globalCoords{ c_passingPoint->getGlobalCoords() };
		arrayS3 globalComponents{ c_normalVector->getLocalComponents() };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			double EC_i = globalComponents[iCoord];
			c_EC[iCoord] = EC_i;
			EC_3 -= EC_i * globalCoords[iCoord];
		}
		c_EC[3] = EC_3;
	}

	/// <summary>
	/// Internal method to apply EC (i.e. create passing point and normal vector)
	/// </summary>
	void Plane::applyEC()
	{
		// Passing point
		double coordX{};
		double coordY{};
		double coordZ{};
		if (std::fabs(c_EC[1]) <= c_toleranceGeneral && std::fabs(c_EC[2]) <= c_toleranceGeneral)
		{
			coordX = -c_EC[3] / c_EC[0];
		}
		else if (std::fabs(c_EC[0]) <= c_toleranceGeneral && std::fabs(c_EC[2]) <= c_toleranceGeneral)
		{
			coordY = -c_EC[3] / c_EC[1];
		}
		else if (std::fabs(c_EC[0]) <= c_toleranceGeneral && std::fabs(c_EC[1]) <= c_toleranceGeneral)
		{
			coordZ = -c_EC[3] / c_EC[2];
		}
		else if (std::fabs(c_EC[0]) <= c_toleranceGeneral)
		{
			coordZ = -c_EC[3] / c_EC[2]; // for y = 0
		}
		else if (std::fabs(c_EC[1]) <= c_toleranceGeneral)
		{
			coordX = -c_EC[3] / c_EC[0]; // for z = 0
		}
		else if (std::fabs(c_EC[2]) <= c_toleranceGeneral)
		{
			coordY = -c_EC[3] / c_EC[1]; // for x = 0
		}
		else
		{
			coordZ = -c_EC[3] / c_EC[2]; // for x = 0 and y = 0
		}
		Point3D passingPoint { arrayS3{ coordX, coordY, coordZ } };
		c_passingPoint = shared_ptr<Point3D>(&passingPoint);

		// The normal vector
		arrayS3 globalComponents;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			globalComponents[iCoord] = c_EC[iCoord];
		}
		Vector3D normalVector { globalComponents };
		c_normalVector = shared_ptr<Vector3D>(&normalVector);
	}
}
