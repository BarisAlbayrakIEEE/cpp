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
	/// The main constructor
	/// </summary>
	/// <exception> ZeroDimensionException: Coincident points </exception>
	Line::Line(
		PointBase& theEndPoint0,
		PointBase& theEndPoint1)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (theEndPoint0.equals(theEndPoint1))
		{
			throw ZeroDimensionException();
		}

		c_endPoint0 = shared_ptr<PointBase>(&theEndPoint0);
		c_endPoint1 = shared_ptr<PointBase>(&theEndPoint1);
		c_length = theEndPoint0.calculateDistance(theEndPoint1);

		// Axis member initialization is not included in the member initializations
		// and performed after invariant inspections (RAII requirement)
		// as it includes object creation which shall be performed after invariant inspections (RAII requirement)
		if (
			theEndPoint0.is2D() &&
			theEndPoint1.is2D() &&
			theEndPoint0.getReferenceCoordSystem() == theEndPoint1.getReferenceCoordSystem()) {
			VectorBase directionVector { DIMENSIONS::D2, Point2D::DownCast(theEndPoint0), Point2D::DownCast(theEndPoint1) };
			Axis axis { theEndPoint0, directionVector };
			c_axis = shared_ptr<Axis>(&axis);
		}
		else {
			VectorBase directionVector { DIMENSIONS::D3, theEndPoint0, theEndPoint1 };
			Axis axis { theEndPoint0, directionVector };
			c_axis = shared_ptr<Axis>(&axis);
		}
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Line::Line(const Line& rhs)
	{
		copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Line& Line::operator=(const Line& rhs) {
		if (&rhs == this) return *this;

		copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Line::Line(Line&& rhs) noexcept
	{
		copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Line& Line::operator=(Line&& rhs) noexcept
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
	bool Line::operator==(const Line& rhs) {
		if (&rhs == this) return true;
		if (*c_endPoint0 != rhs.getEndPoint0()) return false;
		if (*c_endPoint1 != rhs.getEndPoint1()) return false;
		return bool(*c_axis == rhs.getAxis());
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Line::operator!=(const Line& rhs) {
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Line::operator+=(const Line& rhs) {
		if (&rhs == this) return true;
		if (!c_endPoint0->equalsGeometrically(rhs.getEndPoint0())) return false;
		return c_endPoint1->equalsGeometrically(rhs.getEndPoint1());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Line::operator-=(const Line& rhs) {
		return !operator+=(rhs);
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	Line::~Line() {
		Destroy();
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void Line::copyBase(const Line& rhs)
	{
		GeometryObject::copyBase(rhs);
		Axis axis = rhs.getAxis();
		PointBase endPoint0 = rhs.getEndPoint0();
		PointBase endPoint1 = rhs.getEndPoint1();
		
		c_axis = shared_ptr<Axis>(&axis);
		c_endPoint0 = shared_ptr<PointBase>(&endPoint0);
		c_endPoint1 = shared_ptr<PointBase>(&endPoint1);
		c_length = rhs.getLength();
	}

	/// <summary>
	/// Axis, Line, Circle and Plane types do not have additional 2D and 3D types (e.g. Circle2D)
	/// Hence, these types do not have is2D and is3D methods.
	/// They are assumed 3D by default.
	/// See project main docstring in GeometryObject.hxx for more detailed description.
	/// However, these types can be geometrically 2D.
	/// This method defines the conditions to have a 2D object.
	/// </summary>
	bool Line::is2D() const
	{
		return ReferenceObject::is2DStrict(*c_endPoint0, *c_endPoint1);
	}

	/// <summary>
	/// Returns if the instance is not 2D.
	/// See docstring of is2D for details.
	/// </summary>
	bool Line::is3D() const
	{
		return !is2D();
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	void Line::Destroy() {
		c_axis = nullptr;
		c_endPoint0 = nullptr;
		c_endPoint1 = nullptr;
		c_length = 0.;
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool Line::equalsBase(ARGCOPY(Line) theLine) const
	{
		if (this == &theLine)
		{
			return true;
		}
		return true;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool Line::equals(ARGCOPY(Line) theLine) const
	{
		if (!equalsBase(theLine))
		{
			return false;
		}
		if (!c_axis->equals(theLine.getAxis()))
		{
			return false;
		}
		if (!c_endPoint0->equals(theLine.getEndPoint0()))
		{
			return false;
		}
		if (!c_endPoint1->equals(theLine.getEndPoint1()))
		{
			return false;
		}
		return true;
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool Line::equalsGeometrically(ARGCOPY(Line) theLine) const
	{
		if (!equalsBase(theLine))
		{
			return false;
		}
		if (!c_axis->equalsGeometrically(theLine.getAxis()))
		{
			return false;
		}
		if (!c_endPoint0->equalsGeometrically(theLine.getEndPoint0()))
		{
			return false;
		}
		if (!c_endPoint1->equalsGeometrically(theLine.getEndPoint1()))
		{
			return false;
		}
		return true;
	}

	/// <summary>
	/// Getter - Axis
	/// </summary>
	Axis& Line::getAxis() const
	{
		return *c_axis;
	}

	/// <summary>
	/// Getter - Direction vector of the axis
	/// </summary>
	VectorBase& Line::getDirectionVector() const
	{
		return c_axis->getDirectionVector();
	}

	/// <summary>
	/// Getter - End points
	/// </summary>
	std::vector<PointBase> Line::getEndPoints() const
	{
		std::vector<PointBase> outEndPoints;
		outEndPoints.push_back(*c_endPoint0);
		outEndPoints.push_back(*c_endPoint1);
		return outEndPoints;
	}

	/// <summary>
	/// Getter - End point - 0
	/// </summary>
	PointBase& Line::getEndPoint0() const
	{
		return *c_endPoint0;
	}

	/// <summary>
	/// Getter - End point - 1
	/// </summary>
	PointBase& Line::getEndPoint1() const
	{
		return *c_endPoint1;
	}

	/// <summary>
	/// Getter - Length
	/// </summary>
	double Line::getLength() const
	{
		return c_length;
	}

	/// <summary>
	/// Returns the reference CS which is common to all ReferenceObject members.
	///		Returns a handle with nullptr if the ReferenceObject members have different reference CSs.
	/// See module docstring in Axis.hxx for the details.
	/// </summary>
	CoordSystem* Line::getCommonReferenceCoordSystem() const
	{
		if (!is2D())
		{
			return nullptr;
		}
		return &(c_endPoint0->getReferenceCoordSystem());
	}

	/// <summary>
	/// Setter - Direction vector
	/// </summary>
	void Line::setAxis(Axis& theAxis)
	{
		c_axis = shared_ptr<Axis>(&theAxis);
	}

	/// <summary>
	/// Setter - End point 0
	/// </summary>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroDimensionException </exception>
	void Line::setEndPoint0(PointBase& theEndPoint0)
	{
		if (is2D()) {
			if (!GeometryObject::inspectReferenceCoordSystems(theEndPoint0, *c_endPoint1))
			{
				throw CoordSystemMismatchException();
			}
		}
		if (c_endPoint0->equals(*c_endPoint1))
		{
			throw ZeroDimensionException();
		}
		c_endPoint0 = shared_ptr<PointBase>(&theEndPoint0);
		if (c_endPoint1) c_length = c_endPoint0->calculateDistance(*c_endPoint1);
	}

	/// <summary>
	/// Setter - End point1
	/// </summary>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroDimensionException </exception>
	void Line::setEndPoint1(PointBase& theEndPoint1)
	{
		if (is2D()) {
			if (!GeometryObject::inspectReferenceCoordSystems(*c_endPoint0, theEndPoint1))
			{
				throw CoordSystemMismatchException();
			}
		}
		if (c_endPoint1->equals(*c_endPoint0))
		{
			throw ZeroDimensionException();
		}
		c_endPoint1 = shared_ptr<PointBase>(&theEndPoint1);
		if (!c_endPoint0) c_length = c_endPoint0->calculateDistance(*c_endPoint1);
	}

	/// <summary>
	/// Returns if the point is on the line
	/// </summary>
	bool Line::includes(ARGCOPY(PointBase) thePoint) const
	{
		// Distance to axis
		bool isIncluded{ c_axis->includes(thePoint) };
		if (!isIncluded) return false;

		// The two distances from the point to the end points shall be less than the length of the line
		double distance{ thePoint.calculateDistance(*c_endPoint0) };
		if (distance > c_length + c_toleranceGeneral) return false;
		distance = thePoint.calculateDistance(*c_endPoint1);
		if (distance > c_length + c_toleranceGeneral) return false;
		return true;
	}

	/// <summary>
	/// Returns if the line intersects with the input axis
	///		Call the intersects method of the input axis for this line
	/// </summary>
	bool Line::intersects(ARGCOPY(Axis) theAxis) const
	{
		return theAxis.intersects(*this);
	}

	/// <summary>
	/// Returns if the line intersects with the input line
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </summary>
	bool Line::intersects(ARGCOPY(Line) theLine) const
	{
		std::pair<INTERSECTION1, PointBase*> intersectionResults { intersect(theLine) };
		return intersectionResults.first == 0;
	}

	/// <summary>
	/// Returns if this line and the input axis are skew
	///		Call the isSkew method of the input axis for this line
	/// </summary>
	bool Line::isSkew(ARGCOPY(Axis) theAxis) const
	{
		return theAxis.isSkew(*this);
	}

	/// <summary>
	/// Returns if this line and the input line are skew
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </summary>
	bool Line::isSkew(ARGCOPY(Line) theLine) const
	{
		std::pair<INTERSECTION1, PointBase*> intersectionResults { intersect(theLine) };
		return intersectionResults.first == -1;
	}

	/// <summary>
	/// Returns if the line coincides with the input axis
	///		Call the coincides method of the input axis for this line
	/// </sumrnary>
	bool Line::coincides(ARGCOPY(Axis) theAxis) const
	{
		return theAxis.coincides(*this);
	}

	/// <summary>
	/// Returns if the line coincides with the input line
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </sumrnary>
	bool Line::coincides(ARGCOPY(Line) theLine) const
	{
		std::pair<int, PointBase*> intersectionResults{ intersect(theLine) };
		return intersectionResults.first == 1;
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input axis
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	std::pair<INTERSECTION1, PointBase*> Line::intersect(ARGCOPY(Axis) theAxis) const
	{
		return theAxis.intersect(*this);
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	std::pair<INTERSECTION1, PointBase*> Line::intersect(ARGCOPY(Line) theLine) const
	{
		// The axis of the input line
		Axis axis1 { theLine.getAxis() };
		VectorBase& directionVector0 { c_axis->getDirectionVector() };
		VectorBase& directionVector1 { axis1.getDirectionVector() };

		std::pair<INTERSECTION1, PointBase*> intersectionResults{ axis1.intersect(*this) };
		if (intersectionResults.first == INTERSECTION1::Skew1) return intersectionResults;

		// The cross product
		try {
			Vector3D crossProduct = directionVector0.crossProduct(directionVector1);
		} catch(ZeroVectorException) {
			if (includes(theLine.getEndPoint0()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (includes(theLine.getEndPoint1()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (theLine.includes(getEndPoint0()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (theLine.includes(getEndPoint1()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else
			{
				intersectionResults.first = INTERSECTION1::Skew1;
			}
			intersectionResults.second = nullptr;
			return intersectionResults;
		}

		Vector3D crossProduct = directionVector0.crossProduct(directionVector1);
		intersectionResults = c_axis->intersectBase(axis1, crossProduct.getGlobalComponents());
		if (intersectionResults.first != 0)
		{
			return intersectionResults;
		}

		bool isIncluded{ includes(*intersectionResults.second) };
		if (!isIncluded) {
			intersectionResults.first = INTERSECTION1::Skew1;
			return intersectionResults;
		}
		isIncluded = theLine.includes(*intersectionResults.second);
		intersectionResults.first = !isIncluded ? INTERSECTION1::Skew1 : INTERSECTION1::Intersects1;

		return intersectionResults;
	}

	/// <summary>
	/// Projects the input point onto the line
	/// Returns null handle if the projection (onto the axis) is not included by the line
	/// </summary>
	PointBase Line::project(ARGCOPY(PointBase) thePoint)
	{
		if (!includes(thePoint))
		{
			return PointBase(thePoint);
		}
		return PointBase (c_axis->project(thePoint));
	}

	/// <summary>
	/// Calculate distance to a axis
	/// </summary>
	double Line::calculateDistance(ARGCOPY(Axis) theAxis)
	{
		return theAxis.calculateDistance(*this);
	}

	/// <sumrnary>
	/// Calculate distance to a point
	/// </summary>
	double Line::calculateDistance(ARGCOPY(PointBase) thePoint)
	{
		// Get the item in my reference CS
		if (!includes(thePoint)) return 0.;

		PointBase projection{ project(thePoint) };
		double distancePerpandicular{ thePoint.calculateDistance(projection) };
		double distanceProjectionToEndPoint0{ projection.calculateDistance(*c_endPoint0) };
		double distanceProjectionToEndPoint1{ projection.calculateDistance(*c_endPoint1) };
		double distanceProjectionToEndPoint{ std::fmin(distanceProjectionToEndPoint0, distanceProjectionToEndPoint1) };

		return std::pow(std::pow(distancePerpandicular, 2.) + std::pow(distanceProjectionToEndPoint, 2.), 0.5);
	}

	/// <summary>
	/// Calculate distance to a line
	/// </summary>
	double Line::calculateDistance(ARGCOPY(Line) theLine)
	{
		// Inspect if any of the 4 end points is included by the other line
		if (includes(theLine.getEndPoint0())) return 0.;
		if (includes(theLine.getEndPoint1())) return 0.;
		if (theLine.includes(*c_endPoint0)) return 0.;
		if (theLine.includes(*c_endPoint1)) return 0.;

		// Distance to axis
		Axis axis1{ theLine.getAxis() };
		VectorBase& directionVector0{ c_axis->getDirectionVector() };
		VectorBase& directionVector1{ axis1.getDirectionVector() };

		// Distance to axis
		double distanceInfinite{ axis1.calculateDistance(*this) };

		// The axes are parallel
		// Create lines perpandicular to both lines at the end points of the edges (i.e. 4 lines)
		// The distance is same as the axis distance
		// if any of the 4 lines intersects with the other line.
		// Otherwise, the distance is the min of the distances between the endpoints
		if (c_axis->isParallel(axis1)) {
			Vector3D vectorEndToEnd { *c_endPoint0, theLine.getEndPoint0() };
			Vector3D vectorNormal { directionVector0.crossProduct(vectorEndToEnd) };
			Vector3D vectorLineToLine { vectorNormal.crossProduct(directionVector0) };

			Axis axis0 { *c_endPoint0, vectorLineToLine };
			if (axis0.intersects(theLine))
			{
				return distanceInfinite;
			}

			Axis axis1 { *c_endPoint0, vectorLineToLine };
			if (axis1.intersects(theLine))
			{
				return distanceInfinite;
			}

			Axis axis2 { *c_endPoint0, vectorLineToLine };
			if (axis2.intersects(theLine))
			{
				return distanceInfinite;
			}

			Axis axis3 { *c_endPoint0, vectorLineToLine };
			if (axis3.intersects(theLine))
			{
				return distanceInfinite;
			}

			double distance00 = c_endPoint0->calculateDistance(theLine.getEndPoint0());
			double distance01 = c_endPoint0->calculateDistance(theLine.getEndPoint1());
			double distance10 = c_endPoint1->calculateDistance(theLine.getEndPoint0());
			double distance11 = c_endPoint1->calculateDistance(theLine.getEndPoint1());

			return std::fmin(std::fmin(distance00, distance01), std::fmin(distance10, distance11));
		}

		// The axiss intersect
		if (GeometryMath::equals(distanceInfinite, 0., c_toleranceGeneral)) {
			std::pair<int, PointBase*> intersectionResults{ intersect(axis1) };
			if (theLine.includes(*intersectionResults.second)) return 0.;

			std::vector<PointBase> endPoints1{ theLine.getEndPoints() };
			double outDistance{ std::fmin(calculateDistance(endPoints1[0]), calculateDistance(endPoints1[1])) };
			return outDistance;
		}

		// The axiss are skew
		std::vector<PointBase> closestPoints{ c_axis->findClosestPoints(axis1) };
		if (theLine.includes(closestPoints[1]))
		{
			double outDistance{ closestPoints[0].calculateDistance(closestPoints[1]) };
			return outDistance;
		}

		std::vector<PointBase> endPoints1{ theLine.getEndPoints() };
		double outDistance{ std::fmin(
			closestPoints[0].calculateDistance(endPoints1[0]),
			closestPoints[0].calculateDistance(endPoints1[1])) };

		return outDistance;
	}

	/// <summary>
	/// Create a point at the middle of the line
	/// </summary>
	PointBase Line::createMidpoint() const {
		return c_endPoint0->createMidPoint(*c_endPoint1);
	}
}
