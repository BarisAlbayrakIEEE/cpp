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
	/// The main constructor
	/// </summary>
	/// <exception> ZeroDimensionException: Coincident points </exception>
	Line::Line(
		const std::shared_ptr<Point3D>& theEndPoint0,
		const std::shared_ptr<Point3D>& theEndPoint1)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		if (theEndPoint0->equals(*theEndPoint1))
		{
			throw ZeroDimensionException();
		}
		c_endPoint0 = theEndPoint0;
		c_endPoint1 = theEndPoint1;
		c_length = theEndPoint0->calculateDistance(*theEndPoint1);
		c_axis = std::make_shared<Axis>(
			theEndPoint0,
			std::make_shared<Vector3D>(*theEndPoint0, *theEndPoint1));
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool Line::operator==(const Line& rhs) const
	{
		if (&rhs == this) return true;
		if (*c_endPoint0 != *rhs.getEndPoint0()) return false;
		if (*c_endPoint1 != *rhs.getEndPoint1()) return false;
		return bool(*c_axis == *rhs.getAxis());
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool Line::operator!=(const Line& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool Line::operator+=(const Line& rhs) const
	{
		if (&rhs == this) return true;
		if (!c_endPoint0->equalsGeometrically(*rhs.getEndPoint0())) return false;
		return c_endPoint1->equalsGeometrically(*rhs.getEndPoint1());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool Line::operator-=(const Line& rhs) const
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
	/// See == operator docstring
	/// </summary>
	bool Line::equals(ARGCOPY(Line) theLine) const
	{
		if (!c_axis->equals(*theLine.getAxis()))
		{
			return false;
		}
		if (!c_endPoint0->equals(*theLine.getEndPoint0()))
		{
			return false;
		}
		if (!c_endPoint1->equals(*theLine.getEndPoint1()))
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
		if (!c_axis->equalsGeometrically(*theLine.getAxis()))
		{
			return false;
		}
		if (!c_endPoint0->equalsGeometrically(*theLine.getEndPoint0()))
		{
			return false;
		}
		if (!c_endPoint1->equalsGeometrically(*theLine.getEndPoint1()))
		{
			return false;
		}
		return true;
	}

	/// <summary>
	/// Getter - Axis
	/// </summary>
	auto Line::getAxis() const -> std::shared_ptr<Axis>
	{
		return c_axis;
	}

	/// <summary>
	/// Getter - Direction vector of the axis
	/// </summary>
	auto Line::getDirectionVector() const -> std::shared_ptr<Vector3D>
	{
		return c_axis->getDirectionVector();
	}

	/// <summary>
	/// Getter - End point - 0
	/// </summary>
	auto Line::getEndPoint0() const -> std::shared_ptr<Point3D>
	{
		return c_endPoint0;
	}

	/// <summary>
	/// Getter - End point - 1
	/// </summary>
	auto Line::getEndPoint1() const -> std::shared_ptr<Point3D>
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
	auto Line::getCommonReferenceCoordSystem() const -> std::shared_ptr<CoordSystem>
	{
		if (!is2D())
		{
			return nullptr;
		}
		return c_endPoint0->getReferenceCoordSystem();
	}

	/// <summary>
	/// Setter - Direction vector
	/// </summary>
	void Line::setAxis(const std::shared_ptr<Axis>& theAxis)
	{
		c_axis = theAxis;
	}

	/// <summary>
	/// Setter - End point 0
	/// </summary>
	/// <exception> ZeroDimensionException </exception>
	void Line::setEndPoint0(const std::shared_ptr<Point3D>& theEndPoint0)
	{
		if (c_endPoint0->equals(*c_endPoint1))
		{
			throw ZeroDimensionException();
		}
		c_endPoint0 = theEndPoint0;
		if (c_endPoint1) c_length = c_endPoint0->calculateDistance(*c_endPoint1);
	}

	/// <summary>
	/// Setter - End point1
	/// </summary>
	/// <exception> CoordSystemMismatchException </exception>
	/// <exception> ZeroDimensionException </exception>
	void Line::setEndPoint1(const std::shared_ptr<Point3D>& theEndPoint1)
	{
		if (c_endPoint1->equals(*c_endPoint0))
		{
			throw ZeroDimensionException();
		}
		c_endPoint1 = theEndPoint1;
		if (!c_endPoint0) c_length = c_endPoint0->calculateDistance(*c_endPoint1);
	}

	/// <summary>
	/// Returns if the point is on the line
	/// </summary>
	bool Line::includes(ARGCOPY(Point3D) thePoint) const
	{
		// Distance to axis
		if (!c_axis->includes(thePoint)) return false;

		// The two distances from the point to the end points shall be less than the length of the line
		double distance{ thePoint.calculateDistance(*c_endPoint0) };
		if (distance > c_length + getToleranceGeneral()) return false;
		distance = thePoint.calculateDistance(*c_endPoint1);
		if (distance > c_length + getToleranceGeneral()) return false;
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
		return intersect(theLine).first == 0;
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
		return intersect(theLine).first == -1;
	}

	/// <summary>
	/// Returns if the line coincides with the input axis
	///		Call the coincides method of the input axis for this line
	/// </summary>
	bool Line::coincides(ARGCOPY(Axis) theAxis) const
	{
		return theAxis.coincides(*this);
	}

	/// <summary>
	/// Returns if the line coincides with the input line
	///		Call the intersect method which returns the intersection status together with the intersection if exists
	/// </summary>
	bool Line::coincides(ARGCOPY(Line) theLine) const
	{
		return intersect(theLine).first == 1;
	}

	/// <summary>
	/// Returns the intersection status and the intersection point (if exists) with the input axis
	/// Possible cases:
	///		SKew
	///		Intersect
	///		Coincide
	/// </summary>
	auto Line::intersect(ARGCOPY(Axis) theAxis) const
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
	auto Line::intersect(ARGCOPY(Line) theLine) const -> std::pair<INTERSECTION1, std::shared_ptr<Point3D>>
	{
		// The axis of the input line
		auto axis1 { theLine.getAxis() };
		auto directionVector0 { c_axis->getDirectionVector() };
		auto directionVector1 { axis1->getDirectionVector() };

		auto intersectionResults{ axis1->intersect(*this) };
		if (intersectionResults.first == INTERSECTION1::Skew1) return intersectionResults;

		// The cross product
		try {
			auto crossProduct = directionVector0->crossProduct(*directionVector1);
		} catch(ZeroVectorException&) {
			if (includes(*theLine.getEndPoint0()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (includes(*theLine.getEndPoint1()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (theLine.includes(*getEndPoint0()))
			{
				intersectionResults.first = INTERSECTION1::Coincides1;
			}
			else if (theLine.includes(*getEndPoint1()))
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

		auto crossProduct = directionVector0->crossProduct(*directionVector1);
		intersectionResults = c_axis->intersectBase(*axis1, crossProduct->getGlobalComponents());
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
	auto Line::project(ARGCOPY(Point3D) thePoint) const -> std::shared_ptr<Point3D>
	{
		if (includes(thePoint))
		{
			return std::make_shared<Point3D>(thePoint);
		}
		auto project = c_axis->project(thePoint);
		if (!includes(*project.get()))
		{
			return nullptr;
		}
		return project;
	}

	/// <summary>
	/// Calculate distance to a axis
	/// </summary>
	double Line::calculateDistance(ARGCOPY(Axis) theAxis) const
	{
		return theAxis.calculateDistance(*this);
	}

	/// <summary>
	/// Calculate distance to a point
	/// </summary>
	double Line::calculateDistance(ARGCOPY(Point3D) thePoint) const
	{
		// Get the item in my reference CS
		if (!includes(thePoint)) return 0.;

		auto projection{ project(thePoint) };
		double distancePerpandicular{ thePoint.calculateDistance(*projection) };
		double distanceProjectionToEndPoint0{ projection->calculateDistance(*c_endPoint0) };
		double distanceProjectionToEndPoint1{ projection->calculateDistance(*c_endPoint1) };
		double distanceProjectionToEndPoint{ std::fmin(distanceProjectionToEndPoint0, distanceProjectionToEndPoint1) };

		return std::pow(std::pow(distancePerpandicular, 2.) + std::pow(distanceProjectionToEndPoint, 2.), 0.5);
	}

	/// <summary>
	/// Calculate distance to a line
	/// </summary>
	double Line::calculateDistance(ARGCOPY(Line) theLine) const
	{
		// Inspect if any of the 4 end points is included by the other line
		if (includes(*theLine.getEndPoint0())) return 0.;
		if (includes(*theLine.getEndPoint1())) return 0.;
		if (theLine.includes(*c_endPoint0)) return 0.;
		if (theLine.includes(*c_endPoint1)) return 0.;

		// Distance to axis
		auto axis0{ theLine.getAxis() };
		auto directionVector0{ c_axis->getDirectionVector() };
		auto directionVector1{ axis0->getDirectionVector() };

		// Distance to axis
		double distanceInfinite{ calculateDistance(*axis0) };

		// The axes are parallel
		// Create lines perpandicular to both lines at the end points of the edges (i.e. 4 lines)
		// The distance is same as the axis distance
		// if any of the 4 lines intersects with the other line.
		// Otherwise, the distance is the min of the distances between the endpoints
		if (c_axis->isParallel(*axis0)) {
			Vector3D vectorEndToEnd { *c_endPoint0, *theLine.getEndPoint0() };
			auto vectorNormal { directionVector0->crossProduct(vectorEndToEnd) };
			auto vectorLineToLine { vectorNormal->crossProduct(*directionVector0) };

			Axis axis1 { c_endPoint0, vectorLineToLine };
			if (axis1.intersects(theLine))
			{
				return distanceInfinite;
			}

			axis1 = Axis(c_endPoint1, vectorLineToLine);
			if (axis1.intersects(theLine))
			{
				return distanceInfinite;
			}

			axis1 = Axis(theLine.getEndPoint0(), vectorLineToLine);
			if (intersects(axis1))
			{
				return distanceInfinite;
			}

			axis1 = Axis(theLine.getEndPoint1(), vectorLineToLine);
			if (intersects(axis1))
			{
				return distanceInfinite;
			}

			double distance00 = c_endPoint0->calculateDistance(*theLine.getEndPoint0());
			double distance01 = c_endPoint0->calculateDistance(*theLine.getEndPoint1());
			double distance10 = c_endPoint1->calculateDistance(*theLine.getEndPoint0());
			double distance11 = c_endPoint1->calculateDistance(*theLine.getEndPoint1());
			return std::fmin(std::fmin(distance00, distance01), std::fmin(distance10, distance11));
		}

		// The axiss intersect
		if (GeometryMath::equals(distanceInfinite, 0., getToleranceGeneral())) {
			auto intersectionResults{ intersect(*axis0) };
			if (theLine.includes(*intersectionResults.second)) return 0.;

			return std::fmin(
				calculateDistance(*theLine.getEndPoint0()),
				calculateDistance(*theLine.getEndPoint1()));
		}

		// The axiss are skew
		auto closestPoints{ c_axis->findClosestPoints(*axis0) };
		if (theLine.includes(*closestPoints[1]))
		{
			return closestPoints[0]->calculateDistance(*closestPoints[1]);
		}

		return std::fmin(
			closestPoints[0]->calculateDistance(*theLine.getEndPoint0()),
			closestPoints[0]->calculateDistance(*theLine.getEndPoint1()));
	}

	/// <summary>
	/// Create a point at the middle of the line
	/// </summary>
	auto Line::createMidpoint() const -> std::shared_ptr<Point3D>
	{
		if (c_endPoint0->getReferenceCoordSystem()->equals(*c_endPoint1->getReferenceCoordSystem())) {
			auto localCoords = c_endPoint0->createMidPoint(*c_endPoint1)->getLocalCoords();
			return std::make_shared<Point3D>(
					c_endPoint0->getReferenceCoordSystem(),
					localCoords);
		}
		return std::make_shared<Point3D>(c_endPoint0->createMidPoint(*c_endPoint1)->getGlobalCoords());
	}
}
