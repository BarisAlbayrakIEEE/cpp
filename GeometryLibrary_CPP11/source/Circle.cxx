// baris.albayrak.ieee@gmail.com

#include "Macros.h"
#include "GeometryObject.hxx"
#include "GeometryMath.hxx"
#include "GeometryParameters.hxx"
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
	/// The main ctor
	/// Follows RAII idiom
	/// </summary>
	/// <exception> ZeroDimensionException: Radius too small </exception>
	/// <exception> GeometricalMismatchException: The point does not lie on the plane </exception>
	Circle::Circle(
		const std::shared_ptr<Plane>& theReferencePlane,
		const std::shared_ptr<Point3D>& theCenterPoint,
		const double& theRadius)
		:
		GeometryObject()
	{
		if (theRadius < GeometryParameters::getToleranceGeneral())
		{
			throw ZeroDimensionException();
		}
		if (!theReferencePlane->includes(*theCenterPoint))
		{
			throw GeometricalMismatchException();
		}
		c_referencePlane = theReferencePlane;
		c_centerPoint = theCenterPoint;
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
	/// 
	/// Follows RAII idiom
	/// </summary>
	/// <exception> CoincidenceException: At least two points are coincident </exception>
	/// <exception> ColinearPointsException: The input pints lie on the same line </exception>
	/// <exception> UncaughtException: Intersection operation may fail unexpectedly </exception>
	Circle::Circle(
		Point3D thePoint0, // Function needs a copy of the argument
		Point3D thePoint1, // Function needs a copy of the argument
		Point3D thePoint2) // Function needs a copy of the argument
		: GeometryObject()
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

		auto copyPoint0_sp = std::shared_ptr<Point3D>(&thePoint0);
		Axis colinearityLine { copyPoint0_sp, thePoint1 };
		if (colinearityLine.includes(thePoint2))
		{
			throw ColinearPointsException();
		}

		auto copyPoint1_sp = std::shared_ptr<Point3D>(&thePoint1);
		auto copyPoint2_sp = std::shared_ptr<Point3D>(&thePoint2);

		// The reference plane
		Plane referencePlane { copyPoint0_sp, thePoint1, thePoint2 };

		// Step 1
		// Create two lines between the input thePoints
		// (no need to inspect exception as coincidance and colinearity are inspected before)
		Line line0 { copyPoint0_sp, copyPoint1_sp };
		Line line1 { copyPoint1_sp, copyPoint2_sp };

		// Step 2
		// Create thePoints at the middle of the lines
		auto midPoint3D0 = line0.createMidpoint();
		auto midPoint3D1 = line1.createMidpoint();

		Point3D midPoint0 { midPoint3D0->getReferenceCoordSystem(), midPoint3D0->getLocalCoords() };
		Point3D midPoint1 { midPoint3D1->getReferenceCoordSystem(), midPoint3D1->getLocalCoords() };

		// Step 3
		// Create normal vectors to the lines.
		// The vectors lie on the plane containing the lines
		// (no need to inspect exception as thre thePoints define a plane by default)
		std::shared_ptr<Vector3D> vectorNormal { line0.getDirectionVector()->crossProduct(*line1.getDirectionVector()) };
		std::shared_ptr<Vector3D> vector0 { vectorNormal->crossProduct(*line0.getDirectionVector()) };
		std::shared_ptr<Vector3D> vector1 { vectorNormal->crossProduct(*line1.getDirectionVector()) };

		// Step 4
		// Create axes at the line midpoints with the created vectors
		Axis axis0 { std::shared_ptr<Point3D>(&midPoint0), vector0 };
		Axis axis1 { std::shared_ptr<Point3D>(&midPoint1), vector1 };

		// Step 5
		// Create a thePoint at the intersection of the axes
		auto intersectionResults{ axis0.intersect(axis1) };
		if (intersectionResults.first != 0)
		{
			throw UncaughtException();
		}

		// Step 6
		// Calculate the radius and set the members
		double radius{ thePoint0.calculateDistance(*intersectionResults.second) };
		setMembers(std::shared_ptr<Plane>(&referencePlane), intersectionResults.second, radius);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool Circle::operator==(const Circle& rhs) const
	{
		if (&rhs == this) return true;

		double radius1 = c_radius;
		double radius2 = rhs.getRadius();
		if (!GeometryMath::equal_g(radius1, radius2)) return false;
		if (*c_centerPoint != *rhs.getCenterPoint()) return false;
		return bool(*c_referencePlane == *rhs.getReferencePlane());
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Circle::operator!=(const Circle& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Circle::operator+=(const Circle& rhs) const
	{
		if (&rhs == this) return true;

		double radius1 = c_radius;
		double radius2 = rhs.getRadius();
		if (!GeometryMath::equal_g(radius1, radius2)) return false;
		if (*c_centerPoint -= *rhs.getCenterPoint()) return false;
		return *c_referencePlane += *rhs.getReferencePlane();
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Circle::operator-=(const Circle& rhs) const
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
	bool Circle::is2D() const {
		if (
			!c_centerPoint->is2D() &&
			!GeometryMath::zero_g(c_centerPoint->getLocalCoordZ()))
		{
			return false;
		}
		if (c_centerPoint->getReferenceCoordSystem()->getAxisAsVectorZ()->isParallel(*c_referencePlane->getNormalVector()))
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
	/// See += operator docstring
	/// </summary>
	bool Circle::equalsGeometrically(ARGCOPY(Circle) theCircle) const {
		return operator+=(theCircle);
	}

	/// <summary>
	/// Returns null handle if the member is not set yet
	/// </summary>
	auto Circle::getReferencePlane() const -> std::shared_ptr<Plane>
	{
		return c_referencePlane;
	}

	/// <summary>
	/// Returns null handle if the member is not set yet
	/// </summary>
	auto Circle::getCenterPoint() const -> std::shared_ptr<Point3D>
	{
		return c_centerPoint;
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
	auto Circle::getCommonReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>
	{
		auto coordSystem1 = c_centerPoint->getReferenceCoordSystem();
		auto coordSystem2 = c_referencePlane->getCommonReferenceCoordSystem();
		if (!coordSystem2) { return nullptr; }
		if (*coordSystem2 != *coordSystem1)
		{
			return nullptr;
		}
		return coordSystem2;
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary >
	void Circle::setMembers(
		const std::shared_ptr<Plane>& theReferencePlane,
		const std::shared_ptr<Point3D>& theCenterPoint,
		const double& theRadius)
	{
		if (GeometryMath::zero_g(theRadius)) throw ZeroDimensionException();
		c_referencePlane = theReferencePlane;
		c_centerPoint = theCenterPoint;
		c_radius = theRadius;
	}

	/// <summary>
	/// Setter - reference plane
	/// </summary>
	void Circle::setReferencePlane(const std::shared_ptr<Plane>& theReferencePlane)
	{
		if (c_centerPoint && !theReferencePlane->includes(*c_centerPoint))
		{
			throw GeometricalMismatchException();
		}
		c_referencePlane = theReferencePlane;
	}

	/// <summary>
	/// Setter - center point
	/// </summary>
	void Circle::setCenterPoint(const std::shared_ptr<Point3D>& theCenterPoint)
	{
		if (c_referencePlane && !c_referencePlane->includes(*theCenterPoint))
		{
			throw GeometricalMismatchException();
		}
		c_centerPoint = theCenterPoint;
	}

	/// <summary>
	/// Setter - radius
	/// </summary>
	void Circle::setRadius(const double& theRadius)
	{
		if (theRadius < GeometryParameters::getToleranceGeneral()) throw ZeroDimensionException();
		c_radius = theRadius;
	}

	/// <summary>
	/// analyzePoint method inspects if the point is on the circle
	/// Call the method
	/// </summary>
	bool Circle::isPointOn(ARGCOPY(Point3D) thePoint) const {
		if (analyzePoint(thePoint) == 0) return true;
		return false;
	}

	/// <summary>
	/// analyzePoint method inspects if the point is involved by the circle
	/// Call the method
	/// </summary>
	bool Circle::isPointInvolved(ARGCOPY(Point3D) thePoint) const {
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
	int Circle::analyzePoint(ARGCOPY(Point3D) thePoint) const
	{
		// Get the item in my reference CS
		if (!c_referencePlane->includes(thePoint)) return -1;

		double distanceToCenter = c_centerPoint->calculateDistance(thePoint);
		double radius = c_radius;
		if (GeometryMath::equal_g(distanceToCenter, radius)) return 0;
		else if (distanceToCenter < c_radius) return 1;
		return 2;
	}
}
