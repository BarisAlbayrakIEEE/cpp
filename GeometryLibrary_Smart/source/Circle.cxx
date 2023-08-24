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
	/// <exception> ZeroDimensionException: Radius too small </exception>
	/// <exception> GeometricalMismatchException: The point does not lie on the plane </exception>
	Circle::Circle(
		Plane& theReferencePlane,
		PointBase& theCenterPoint,
		const double& theRadius)
		:
		GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (theRadius < c_toleranceGeneral)
		{
			throw ZeroDimensionException();
		}
		if (!theReferencePlane.includes(theCenterPoint))
		{
			throw GeometricalMismatchException();
		}
		c_referencePlane = std::shared_ptr<Plane>(&theReferencePlane);
		c_centerPoint = std::shared_ptr<PointBase>(&theCenterPoint);
		c_radius = theRadius;
	}

	/// <summary>
	/// Preconditions:
	///     Points shall not be coincident
	///     Points shall not be co-linear

	/// Theory:
	/// Consider a point pair PP on the circle.
	/// Consider a line (Line A) passing through the points of PP.
	/// Consider an axis (Axis B) which is normal to Line A at its midpoint.
	/// Axis B passes through the circle's centroid for every point pair on the circle
	/// as the length of any line from the centroid to a point on the circle is the circle radius
	/// which means a triangle formed by the circle center and PP is equilateral
	/// which means the height is the median at the same time.

	/// Approach:
	/// Step 1:
	/// Create two lines between the input points
	/// Step 2:
	/// Create points at the middle of the lines
	/// Step 3:
	/// Create normal vectors to the lines. The vectors lie on the plane containing the lines
	/// Step 4:
	/// Create axes at the line midpoints with the created vectors
	/// Step 5:
	/// Create a point at the intersection of the axiss
	/// Step 6:
	/// Calculate the radius and set the members
	/// </summary>
	/// <exception> CoincidenceException: At least two points are coincident </exception>
	/// <exception> ColinearPointsException: The input pints lie on the same line </exception>
	/// <exception> UncaughtException: Intersection operation may fail unexpectedly </exception>
	Circle::Circle(
		PointBase& thePoint0,
		PointBase& thePoint1,
		PointBase& thePoint2)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		// Inspect inputs
		if (thePoint0.coincides(thePoint1))
		{
			throw CoincidenceException();
		}
		if (thePoint0.coincides(thePoint2))
		{
			throw CoincidenceException();
		}

		Axis colinearityLine { thePoint0, thePoint1 };
		if (colinearityLine.includes(thePoint2))
		{
			throw ColinearPointsException();
		}

		// The reference plane
		Plane referencePlane { thePoint0, thePoint1, thePoint2 };

		// Step 1
		// Create two lines between the input thePoints
		// (no need to inspect exception as coincidance and colinearity are inspected before)
		Line line0 { thePoint0, thePoint1 };
		Line line1 { thePoint1, thePoint2 };

		// Step 2
		// Create thePoints at the middle of the lines
		PointBase midPointBase0 = line0.createMidpoint();
		PointBase midPointBase1 = line1.createMidpoint();

		Point3D midPoint0 { midPointBase0.getReferenceCoordSystem(), midPointBase0.getLocalCoords() };
		Point3D midPoint1 { midPointBase1.getReferenceCoordSystem(), midPointBase1.getLocalCoords() };

		// Step 3
		// Create normal vectors to the lines.
		// The vectors lie on the plane containing the lines
		// (no need to inspect exception as thre thePoints define a plane by default)
		Vector3D vectorNormal { line0.getDirectionVector().crossProduct(line1.getDirectionVector()) };
		Vector3D vector0 { vectorNormal.crossProduct(line0.getDirectionVector()) };
		Vector3D vector1 { vectorNormal.crossProduct(line1.getDirectionVector()) };

		// Step 4
		// Create axes at the line midpoints with the created vectors
		Axis axis0 { midPoint0, vector0 };
		Axis axis1 { midPoint1, vector1 };

		// Step 5
		// Create a thePoint at the intersection of the axiss
		std::pair<int, PointBase*> intersectionResults{ axis0.intersect(axis1) };
		if (intersectionResults.first != 0)
		{
			throw UncaughtException();
		}
		Point3D centerPoint = *(static_cast<Point3D*>(intersectionResults.second));

		// Step 6
		// Calculate the radius and set the members
		double radius{ thePoint0.calculateDistance(centerPoint) };
		setMembers(referencePlane, centerPoint, radius);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Circle::Circle(const Circle& rhs)
	{
		copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Circle& Circle::operator=(const Circle& rhs) {
		if (&rhs == this) return *this;

		copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Circle::Circle(Circle&& rhs) noexcept
	{
		copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	Circle& Circle::operator=(Circle&& rhs) noexcept {
		if (&rhs == this) return *this;

		copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool Circle::operator==(const Circle& rhs) {
		if (&rhs == this) return true;
		if (!GeometryMath::equals(c_radius, rhs.getRadius(), c_toleranceGeneral)) return false;
		if (*c_centerPoint != rhs.getCenterPoint()) return false;
		return bool(*c_referencePlane == rhs.getReferencePlane());
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Circle::operator!=(const Circle& rhs) {
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Circle::operator+=(const Circle& rhs) {
		if (&rhs == this) return true;
		if (!GeometryMath::equals(c_radius, rhs.getRadius(), c_toleranceGeneral)) return false;
		if (!c_centerPoint->equalsGeometrically(rhs.getCenterPoint())) return false;
		return c_referencePlane->equalsGeometrically(rhs.getReferencePlane());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Circle::operator-=(const Circle& rhs) {
		return !operator+=(rhs);
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	Circle::~Circle() {
		Destroy();
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void Circle::copyBase(const Circle& rhs)
	{
		GeometryObject::copyBase(rhs);
		Plane referencePlane { rhs.getReferencePlane() };
		PointBase centerPoint { rhs.getCenterPoint() };
		c_referencePlane = std::shared_ptr<Plane>(&referencePlane);
		c_centerPoint = std::shared_ptr<PointBase>(&centerPoint);
		c_radius = rhs.getRadius();
	}

	/// <summary>
	/// Axis, Line, Circle and Plane types do not have additional 2D and 3D types (e.g. Circle2D)
	/// Hence, these types do not have is2D and is3D methods.
	/// They are assumed 3D by default.
	/// See project main docstring in GeometryObject.hxx for more detailed description.
	/// However, these types can be geometrically 2D.
	/// This method defines the conditions to have a 2D object.
	/// </summary>
	bool Circle::is2D() const {
		if (
			!c_centerPoint->is2D() &&
			!GeometryMath::equals(c_centerPoint->getLocalCoordZ(), 0., c_toleranceGeneral))
		{
			return false;
		}
		if (c_centerPoint->getReferenceCoordSystem().getAxisAsVectorZ().isParallel(c_referencePlane->getNormalVector()))
		{
			return true;
		}
		return false;
	}

	/// <summary>
	/// Returns if the instance is not 2D.
	/// See docstring of is2D for details.
	/// </summary>
	bool Circle::is3D() const
	{
		return !is2D();
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	void Circle::Destroy() {
		c_referencePlane = nullptr;
		c_centerPoint = nullptr;
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool Circle::equalsBase(ARGCOPY(Circle) theCircle) const {
		return true;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool Circle::equals(ARGCOPY(Circle) theCircle) const {
		if (this == &theCircle) return true;
		if (!equalsBase(theCircle)) return false;
		if (!c_referencePlane->equals(theCircle.getReferencePlane())) return false;
		if (!c_centerPoint->equals(theCircle.getCenterPoint())) return false;
		return GeometryMath::equals(c_radius, theCircle.getRadius(), c_toleranceGeneral);
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool Circle::equalsGeometrically(ARGCOPY(Circle) theCircle) const {
		if (this == &theCircle) return true;
		if (!equalsBase(theCircle)) return false;
		if (!c_referencePlane->equals(theCircle.getReferencePlane())) return false;
		if (!c_centerPoint->equalsGeometrically(theCircle.getCenterPoint())) return false;
		return GeometryMath::equals(c_radius, theCircle.getRadius(), c_toleranceGeneral);
	}

	/// <summary>
	/// Returns null handle if the member is not set yet
	/// </summary>
	Plane& Circle::getReferencePlane() const {
		return *c_referencePlane;
	}

	/// <summary>
	/// Returns null handle if the member is not set yet
	/// </summary>
	PointBase& Circle::getCenterPoint() const {
		return *c_centerPoint;
	}

	/// </summary>
	/// Returns -1 if the member is not set yet
	/// </summary>
	double Circle::getRadius() const {
		return c_radius;
	}

	/// <summary>
	/// Returns the reference CS which is common to all ReferenceObject members.
	///		Returns a handle with nullptr if the ReferenceObject members have different reference CSs.
	/// See module docstring in Axis.hxx for the details.
	/// </summary>
	/// <exception> ZeroVectorException </exception>
	CoordSystem* Circle::getCommonReferenceCoordSystem() const
	{
		CoordSystem coordSystem1 = c_centerPoint->getReferenceCoordSystem();
		CoordSystem* coordSystem2 = c_referencePlane->getCommonReferenceCoordSystem();
		if (!coordSystem2) { return nullptr; }
		if (!coordSystem2->equals(coordSystem1))
		{
			return nullptr;
		}
		return coordSystem2;
	}

	/// <sumrnary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </sumrnary >
	void Circle::setMembers(
		Plane& theReferencePlane,
		PointBase& theCenterPoint,
		const double& theRadius)
	{
		if (theRadius < c_toleranceGeneral) throw ZeroDimensionException();
		c_referencePlane = std::shared_ptr<Plane>(&theReferencePlane);
		c_centerPoint = std::shared_ptr<PointBase>(&theCenterPoint);
		c_radius = theRadius;
	}

	/// <summary>
	/// Setter - reference plane
	/// </summary>
	void Circle::setReferencePlane(Plane& theReferencePlane)
	{
		if (c_centerPoint) {
			if (!theReferencePlane.includes(*c_centerPoint))
			{
				throw GeometricalMismatchException();
			}
		}
		c_referencePlane = std::shared_ptr<Plane>(&theReferencePlane);
	}

	/// <summary>
	/// Setter - center point
	/// </summary>
	void Circle::setCenterPoint(PointBase& theCenterPoint)
	{
		if (c_referencePlane) {
			if (!c_referencePlane->includes(theCenterPoint))
			{
				throw GeometricalMismatchException();
			}
		}
		c_centerPoint = std::shared_ptr<PointBase>(&theCenterPoint);
	}

	/// <summary>
	/// Setter - radius
	/// </summary>
	void Circle::setRadius(const double& theRadius)
	{
		if (theRadius < c_toleranceGeneral) throw ZeroDimensionException();
		c_radius = theRadius;
	}

	/// <summary>
	/// analyzePoint method inspects if the point is on the circle
	/// Call the method
	/// </summary>
	bool Circle::isPointOn(ARGCOPY(PointBase) thePoint) const {
		if (analyzePoint(thePoint) == 0) return true;
		return false;
	}

	/// <summary>
	/// analyzePoint method inspects if the point is involved by the circle
	/// Call the method
	/// </summary>
	bool Circle::isPointInvolved(ARGCOPY(PointBase) thePoint) const {
		if (analyzePoint(thePoint) == 1) return true;
		return false;
	}

	/// <summary>
	/// Followings are the possible outputs of this method
	///	-1: Point is not on the plane of the circle
	///	0: Point is on the circle
	///	1: Point is on the plane of the circle and involved by the circle
	///	2: Point is on the plane olthe circle but aut olthe circle area
	/// </summary>
	int Circle::analyzePoint(ARGCOPY(PointBase) thePoint) const
	{
		// Get the item in my reference CS
		if (!c_referencePlane->includes(thePoint)) return -1;

		double distanceToCenter = c_centerPoint->calculateDistance(thePoint);
		if (GeometryMath::equals(distanceToCenter, c_radius, c_toleranceGeneral)) return 0;
		else if (distanceToCenter < c_radius) return 1;
		else return 2;
	}
}
