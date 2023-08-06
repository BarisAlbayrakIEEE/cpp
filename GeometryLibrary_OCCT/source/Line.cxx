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
	IMPLEMENT_STANDARD_RTTIEXT(Line, GeometryObject)

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> ZeroDimensionException: Coincident points </exception>
	Line::Line(
		ARGCOPY(PointBase) theEndPoint0,
		ARGCOPY(PointBase) theEndPoint1) throw(NullptrException, ZeroDimensionException)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (theEndPoint0.IsNull())
		{
			throw NullptrException();
		}
		if (theEndPoint1.IsNull())
		{
			throw NullptrException();
		}
		if (*theEndPoint0.get() += *theEndPoint1.get())
		{
			throw ZeroDimensionException();
		}

		c_endPoint0 = theEndPoint0;
		c_endPoint1 = theEndPoint1;
		c_length = c_endPoint0->calculateDistance(c_endPoint1);

		// Axis member initialization is not included in the member initializations
		// and performed after invariant inspections (RAII requirement)
		// as it includes object creation which shall be performed after invariant inspections (RAII requirement)
		Handle(VectorBase) directionVector;
		if (
			theEndPoint0->is2D() &&
			theEndPoint1->is2D() &&
			theEndPoint0->getReferenceCoordSystem() == theEndPoint1->getReferenceCoordSystem()) {
			directionVector = Handle(Vector2D)(
				new Vector2D(
					Handle(Point2D)::DownCast(theEndPoint0), 
					Handle(Point2D)::DownCast(theEndPoint1)));
		}
		else {
			directionVector = Handle(Vector3D)(new Vector3D(theEndPoint0, theEndPoint1));
		}
		c_axis = new Axis(theEndPoint0, directionVector);
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
		if (c_endPoint0 != rhs.getEndPoint0()) return false;
		if (c_endPoint1 != rhs.getEndPoint1()) return false;
		return bool(c_axis == rhs.getAxis());
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
	/// Use OCCT approach
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
		c_axis = rhs.getAxis();
		c_endPoint0 = rhs.getEndPoint0();
		c_endPoint1 = rhs.getEndPoint1();
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
		return ReferenceObject::is2DStrict(c_endPoint0, c_endPoint1);
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
	/// Use Nullify method of the OCCT Standard_Handle for the object destruction
	/// </summary>
	void Line::Destroy() {
		c_axis.Nullify();
		c_endPoint0.Nullify();
		c_endPoint1.Nullify();
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool Line::equalsBase(ARGCOPY(Line) theLine) const
	{
		if (this == theLine.get())
		{
			return true;
		}
		if (theLine.IsNull())
		{
			return false;
		}
		if (!inspectNullEquality(c_axis, theLine->getAxis()))
		{
			return false;
		}
		if (!inspectNullEquality(c_endPoint0, theLine->getEndPoint0()))
		{
			return false;
		}
		if (!inspectNullEquality(c_endPoint1, theLine->getEndPoint1()))
		{
			return false;
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
		if (!c_axis->equals(theLine->getAxis()))
		{
			return false;
		}
		if (!c_endPoint0->equals(theLine->getEndPoint0()))
		{
			return false;
		}
		if (!c_endPoint1->equals(theLine->getEndPoint1()))
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
		if (!c_axis->equalsGeometrically(theLine->getAxis()))
		{
			return false;
		}
		if (!c_endPoint0->equalsGeometrically(theLine->getEndPoint0()))
		{
			return false;
		}
		if (!c_endPoint1->equalsGeometrically(theLine->getEndPoint1()))
		{
			return false;
		}
		return true;
	}

	/// <summary>
	/// Getter - Axis
	/// </summary>
	OUTVAL(Axis) Line::getAxis() const
	{
		return c_axis;
	}

	/// <summary>
	/// Getter - Direction vector of the axis
	/// </summary>
	OUTVAL(VectorBase) Line::getDirectionVector() const
	{
		return c_axis->getDirectionVector();
	}

	/// <summary>
	/// Getter - End points
	/// </summary>
	std::vector<OUTVAL(PointBase)> Line::getEndPoints() const
	{
		std::vector<Handle(PointBase)> outEndPoints;
		outEndPoints.push_back(c_endPoint0);
		outEndPoints.push_back(c_endPoint1);
		return outEndPoints;
	}

	/// <summary>
	/// Getter - End point - 0
	/// </summary>
	OUTVAL(PointBase) Line::getEndPoint0() const
	{
		return c_endPoint0;
	}

	/// <summary>
	/// Getter - End point - 1
	/// </summary>
	OUTVAL(PointBase) Line::getEndPoint1() const
	{
		return c_endPoint1;
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
	/// <exception> ZeroVectorException </exception>
	OUTVAL(CoordSystem) Line::getCommonReferenceCoordSystem() const
	{
		if (!is2D())
		{
			return Handle(CoordSystem)(0);
		}
		return c_endPoint0->getReferenceCoordSystem();
	}

	/// <summary>
	/// Setter - Direction vector
	/// </summary>
	/// <exception> NullptrException </exception>
	void Line::setAxis(ARGCOPY(Axis) theAxis) throw(NullptrException)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		c_axis = theAxis;
	}

	/// <summary>
	/// Setter - End point 0
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroDimensionException </exception>
	void Line::setEndPoint0(ARGCOPY(PointBase) theEndPoint0) throw(
		NullptrException,
		CoordSystemMismatchException,
		ZeroDimensionException)
	{
		if (theEndPoint0.IsNull())
		{
			throw NullptrException();
		}
		if (is2D()) {
			if (!GeometryObject::inspectReferenceCoordSystems(theEndPoint0, c_endPoint1))
			{
				throw CoordSystemMismatchException();
			}
		}
		if (c_endPoint0->equals(c_endPoint1))
		{
			throw ZeroDimensionException();
		}
		c_endPoint0 = theEndPoint0;
		if (!c_endPoint1.IsNull()) c_length = c_endPoint0->calculateDistance(c_endPoint1);
	}

	/// <summary>
	/// Setter - End point1
	/// </summary>
	/// <exception> NullptrException </exception>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroDimensionException </exception>
	void Line::setEndPoint1(ARGCOPY(PointBase) theEndPoint1) throw(
		NullptrException,
		CoordSystemMismatchException,
		ZeroDimensionException)
	{
		if (theEndPoint1.IsNull())
		{
			throw NullptrException();
		}
		if (is2D()) {
			if (!GeometryObject::inspectReferenceCoordSystems(c_endPoint0, theEndPoint1))
			{
				throw CoordSystemMismatchException();
			}
		}
		if (c_endPoint1->equals(c_endPoint0))
		{
			throw ZeroDimensionException();
		}
		c_endPoint1 = theEndPoint1;
		if (!c_endPoint0.IsNull()) c_length = c_endPoint0->calculateDistance(c_endPoint1);
	}

	/// <summary>
	/// Returns if the point is on the line
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Line::includes(ARGCOPY(PointBase) thePoint) throw(NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Distance to axis
		bool isIncluded{ c_axis->includes(thePoint) };
		if (!isIncluded) return false;

		// The two distances from the point to the end points shall be less than the length of the line
		double distance{ thePoint->calculateDistance(c_endPoint0) };
		if (distance > c_length + c_toleranceGeneral) return false;
		distance = thePoint->calculateDistance(c_endPoint1);
		if (distance > c_length + c_toleranceGeneral) return false;
		return true;
	}

	/// <summary>
	/// Returns if the line intersects with the input axis
	///		Call the intersects method of the input axis for this line
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Line::intersects(ARGCOPY(Axis) theAxis) throw(NullptrException)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return theAxis->intersects(this);
	}

	/// <summary>
	/// Returns if the line intersects with the input line
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Line::intersects(ARGCOPY(Line) theLine) throw(NullptrException)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}
		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theLine) };
		return intersectionResults.first == 0;
	}

	/// <summary>
	/// Returns if this line and the input axis are skew
	///		Call the isSkew method of the input axis for this line
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Line::isSkew(ARGCOPY(Axis) theAxis) throw(NullptrException)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return theAxis->isSkew(this);
	}

	/// <summary>
	/// Returns if this line and the input line are skew
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </summary>
	/// <exception> NullptrException </exception>
	bool Line::isSkew(ARGCOPY(Line) theLine) throw(NullptrException)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}
		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theLine) };
		return intersectionResults.first == -1;
	}

	/// <summary>
	/// Returns if the line coincides with the input axis
	///		Call the coincides method of the input axis for this line
	/// </sumrnary>
	/// <exception> NullptrException </exception>
	bool Line::coincides(ARGCOPY(Axis) theAxis) throw(NullptrException)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return theAxis->coincides(this);
	}

	/// <summary>
	/// Returns if the line coincides with the input line
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </sumrnary>
	/// <exception> NullptrException </exception>
	bool Line::coincides(ARGCOPY(Line) theLine) throw(NullptrException)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}
		std::pair<int, Handle(PointBase)> intersectionResults{ intersect(theLine) };
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
	std::pair<INTERSECTION1, Handle(PointBase)> Line::intersect(ARGCOPY(Axis) theAxis)
		throw(NullptrException)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return theAxis->intersect(this);
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input line
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	/// <exception> NullptrException </exception>
	std::pair<INTERSECTION1, Handle(PointBase)> Line::intersect(ARGCOPY(Line) theLine)
		throw(NullptrException)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		// The axis of the input line
		Handle(Axis) axis1{ theLine->getAxis() };
		Handle(VectorBase) directionVector0{ c_axis->getDirectionVector() };
		Handle(VectorBase) directionVector1{ axis1->getDirectionVector() };

		std::pair<INTERSECTION1, Handle(PointBase)> intersectionResults{ axis1->intersect(this) };
		if (intersectionResults.first == INTERSECTION1::Skew1) return intersectionResults;

		// The cross product
		bool chkParallelism{};
		Handle(Vector3D) crossProduct = directionVector0->crossProduct(directionVector1);
		if (crossProduct.IsNull()) {
			chkParallelism = true;
		}
		if (chkParallelism || crossProduct.IsNull() || GeometryMath::equals(crossProduct->getMagnitude(), 0., c_toleranceGeneral)) {
			if (includes(theLine->getEndPoint0()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (includes(theLine->getEndPoint1()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (theLine->includes(getEndPoint0()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (theLine->includes(getEndPoint1()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else
			{
				intersectionResults.first = INTERSECTION1::Skew1;
			}
			return intersectionResults;
		}

		intersectionResults = c_axis->intersectBase(axis1, crossProduct->getGlobalComponents());
		if (intersectionResults.first != 0)
		{
			return intersectionResults;
		}

		bool isIncluded{ includes(intersectionResults.second) };
		if (!isIncluded) {
			intersectionResults.first = INTERSECTION1::Skew1;
			return intersectionResults;
		}
		isIncluded = theLine->includes(intersectionResults.second);
		intersectionResults.first = !isIncluded ? INTERSECTION1::Skew1 : INTERSECTION1::Intersects1;

		return intersectionResults;
	}

	/// <summary>
	/// Projects the input point onto the line
	/// Returns null handle if the projection (onto the axis) is not included by the line
	/// </summary>
	/// <exception> NullptrException </exception>
	OUTVAL(PointBase) Line::project(ARGCOPY(PointBase) thePoint) throw(NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		Handle(PointBase) projection{ c_axis->project(thePoint) };
		if (!includes(projection))
		{
			return Handle(PointBase)(0);
		}
		return projection;
	}

	/// <summary>
	/// Calculate distance to a axis
	/// </summary>
	/// <exception> NullptrException </exception>
	double Line::calculateDistance(ARGCOPY(Axis) theAxis) throw(NullptrException)
	{
		if (theAxis.IsNull())
		{
			throw NullptrException();
		}
		return theAxis->calculateDistance(this);
	}

	/// <sumrnary>
	/// Calculate distance to a point
	/// </summary>
	/// <exception> NullptrException </exception>
	double Line::calculateDistance(ARGCOPY(PointBase) thePoint) throw(NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Get the item in my reference CS
		Handle(PointBase) projection{ project(thePoint) };
		if (projection.IsNull()) return thePoint->calculateDistance(projection);

		projection = getAxis()->project(thePoint);
		double distancePerpandicular{ thePoint->calculateDistance(projection) };
		double distanceProjectionToEndPoint0{ projection->calculateDistance(c_endPoint0) };
		double distanceProjectionToEndPoint1{ projection->calculateDistance(c_endPoint1) };
		double distanceProjectionToEndPoint{ std::fmin(distanceProjectionToEndPoint0, distanceProjectionToEndPoint1) };

		return std::pow(std::pow(distancePerpandicular, 2.) + std::pow(distanceProjectionToEndPoint, 2.), 0.5);
	}

	/// <summary>
	/// Calculate distance to a line
	/// </summary>
	/// <exception> NullptrException </exception>
	double Line::calculateDistance(ARGCOPY(Line) theLine) throw(NullptrException)
	{
		if (theLine.IsNull())
		{
			throw NullptrException();
		}

		// Inspect if any of the 4 end points is included by the other line
		if (includes(theLine->getEndPoint0())) return 0.;
		if (includes(theLine->getEndPoint1())) return 0.;
		if (theLine->includes(c_endPoint0)) return 0.;
		if (theLine->includes(c_endPoint1)) return 0.;

		// Distance to axis
		Handle(Axis) axis1{ theLine->getAxis() };
		Handle(VectorBase) directionVector0{ c_axis->getDirectionVector() };
		Handle(VectorBase) directionVector1{ axis1->getDirectionVector() };

		// Distance to axis
		double distanceInfinite{ axis1->calculateDistance(this) };

		// The axes are parallel
		// Create lines perpandicular to both lines at the end points of the edges (i.e. 4 lines)
		// The distance is same as the axis distance
		// if any of the 4 lines intersects with the other line.
		// Otherwise, the distance is the min of the distances between the endpoints
		if (c_axis->isParallel(axis1)) {
			Handle(Vector3D) vectorEndToEnd { new Vector3D(c_endPoint0, theLine->getEndPoint0()) };
			Handle(Vector3D) vectorNormal { directionVector0->crossProduct(vectorEndToEnd) };
			Handle(Vector3D) vectorLineToLine { vectorNormal->crossProduct(directionVector0) };

			Handle(Axis) axis0{ new Axis(c_endPoint0, vectorLineToLine) };
			if (axis0->intersects(theLine))
			{
				return distanceInfinite;
			}

			Handle(Axis) axis1 { new Axis(c_endPoint0, vectorLineToLine) };
			if (axis1->intersects(theLine))
			{
				return distanceInfinite;
			}

			Handle(Axis) axis2 { new Axis(c_endPoint0, vectorLineToLine) };
			if (axis2->intersects(theLine))
			{
				return distanceInfinite;
			}

			Handle(Axis) axis3 { new Axis(c_endPoint0, vectorLineToLine) };
			if (axis3->intersects(theLine))
			{
				return distanceInfinite;
			}

			double distance00 = c_endPoint0->calculateDistance(theLine->getEndPoint0());
			double distance01 = c_endPoint0->calculateDistance(theLine->getEndPoint1());
			double distance10 = c_endPoint1->calculateDistance(theLine->getEndPoint0());
			double distance11 = c_endPoint1->calculateDistance(theLine->getEndPoint1());

			return std::fmin(std::fmin(distance00, distance01), std::fmin(distance10, distance11));
		}

		// The axiss intersect
		if (GeometryMath::equals(distanceInfinite, 0., c_toleranceGeneral)) {
			std::pair<int, Handle(PointBase)> intersectionResults{ intersect(axis1) };
			if (theLine->includes(intersectionResults.second)) return 0.;

			std::vector<Handle(PointBase)> endPoints1{ theLine->getEndPoints() };
			double outDistance{ std::fmin(calculateDistance(endPoints1[0]), calculateDistance(endPoints1[1])) };
			return outDistance;
		}

		// The axiss are skew
		std::vector<Handle(PointBase)> closestPoints{ c_axis->findClosestPoints(axis1) };
		if (theLine->includes(closestPoints[1]))
		{
			double outDistance{ closestPoints[0]->calculateDistance(closestPoints[1]) };
			return outDistance;
		}

		std::vector<Handle(PointBase)> endPoints1{ theLine->getEndPoints() };
		double outDistance{ std::fmin(
			closestPoints[0]->calculateDistance(endPoints1[0]),
			closestPoints[0]->calculateDistance(endPoints1[1])) };

		return outDistance;
	}

	/// <summary>
	/// Create a point at the middle of the line
	/// </summary>
	OUTVAL(PointBase) Line::createMidpoint() const {
		return c_endPoint0->createMidPoint(c_endPoint1);
	}
}
