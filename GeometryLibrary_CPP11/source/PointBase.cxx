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
	/// The default constructor
	/// Reference CS is the global CS
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(const int theDimensionCount)
		:
		ReferenceObject(theDimensionCount) { }

	/// <summary>
	/// Ctor
	/// Point by coords
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		:
		ReferenceObject(theDimensionCount, theReferenceCoordSystem) { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		const std::array<double, 3>& theLocalCoords)
		:
		ReferenceObject(theDimensionCount),
		c_localCoords{ theLocalCoords } { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		const std::vector<double, std::allocator<double>>& theLocalCoords)
		:
		ReferenceObject(theDimensionCount)
	{
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalCoords)
		:
		ReferenceObject(theDimensionCount, theReferenceCoordSystem),
		c_localCoords{ theLocalCoords } { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalCoords)
		:
		ReferenceObject(theDimensionCount, theReferenceCoordSystem)
	{
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool PointBase::operator==(const PointBase& rhs) const
	{
		if (&rhs == this) return true;
		if (ReferenceObject::operator!=(rhs)) return false;
		return GeometryMath::equals(c_localCoords, rhs.getLocalCoords(), getToleranceGeneral());
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool PointBase::operator!=(const PointBase& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool PointBase::operator+=(const PointBase& rhs) const
	{
		if (&rhs == this) return true;
		return GeometryMath::equals(getGlobalCoords(), rhs.getGlobalCoords(), getToleranceGeneral());
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool PointBase::operator-=(const PointBase& rhs) const
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool PointBase::equals(ARGCOPY(PointBase) thePoint) const {
		if (this == &thePoint) return true;
		if (!ReferenceObject::equalsRef(thePoint)) return false;
		return GeometryMath::equals(c_localCoords, thePoint.getLocalCoords(), getToleranceGeneral());
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool PointBase::equalsGeometrically(ARGCOPY(PointBase) thePoint) const {
		if (this == &thePoint) return true;
		return GeometryMath::equals(
			c_localCoords,
			getPointWithMyCoordSystem(thePoint)->getLocalCoords(),
			getToleranceGeneral());
	}

	/// <summary>
	/// Getter - Local coord - X
	/// </summary>
	double PointBase::getLocalCoordX() const {
		return c_localCoords[0];
	}

	/// <summary>
	/// Getter - Local coord - Y
	/// </summary>
	double PointBase::getLocalCoordY() const {
		return c_localCoords[1];
	}

	/// <summary>
	/// Getter - Local coord - Z
	/// </summary>
	double PointBase::getLocalCoordZ() const {
		return c_localCoords[2];
	}

	/// <summary>
	/// Getter - Local coords
	/// </summary>
	auto PointBase::getLocalCoords() const -> std::array<double, 3>
	{
		return c_localCoords;
	}

	/// <summary>
	/// Getter - Global coord - X
	/// </summary>
	double PointBase::getGlobalCoordX() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localCoords[0];
		return getGlobalCoords()[0];
	}

	/// <summary>
	/// Getter - Global coord - Y
	/// </summary>
	double PointBase::getGlobalCoordY() const {
		if (c_referenceCoordSystem->isGlobal()) return c_localCoords[1];
		return getGlobalCoords()[1];
	}

	/// <summary>
	/// Getter - Global coord - Z
	/// </summary>
	double PointBase::getGlobalCoordZ() const
	{
		if (c_referenceCoordSystem->isGlobal()) return c_localCoords[2];
		return getGlobalCoords()[2];
	}

	/// <summary>
	/// Getter - Global coords
	/// Returns null handle if the localcoords member is not set yet
	/// The global coordinates are the components of the vector:
	///	The vector from global origin to the origin of the reference CS
	///	+
	///	the vector from the origin of the reference CS to the point
	/// </summary>
	auto PointBase::getGlobalCoords() const -> std::array<double, 3>
	{
		if (c_referenceCoordSystem->isGlobal()) return c_localCoords;

		auto axesVectors{ c_referenceCoordSystem->getAxesAsVector() }; // By definition, the reference CS is the global CS
		auto outGlobalCoords{ c_referenceCoordSystem->getOriginCoords() };
		for (int iAxis = 0; iAxis < DIMENSIONS::D3; iAxis++) {
			auto vectorComponents{ axesVectors[iAxis]->getLocalComponents() };
			for (int iComponent = 0; iComponent < DIMENSIONS::D3; iComponent++)
			{
				outGlobalCoords[iComponent] += vectorComponents[iComponent] * c_localCoords[iComponent];
			}
		}
		return outGlobalCoords;
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void PointBase::setMembers(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::array<double, 3>& theLocalCoords)
	{
		ReferenceObject::setMembersRef(theDimensionCount, theReferenceCoordSystem);
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void PointBase::setMembers(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		const std::vector<double, std::allocator<double>>& theLocalCoords)
	{
		ReferenceObject::setMembersRef(theDimensionCount, theReferenceCoordSystem);
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// Setter - reference CS
	/// Public method for the users.
	/// theKeepGlobalCoordsSame:
	///		true:
	///			Keep the coords wrt the GLOBAL CS same and update the coords wrt the REFERENCE CS.
	///		false:
	///			Keep the coords wrt the REFERENCE CS same and update the coords wrt the GLOBAL CS.
	/// </summary>
	void PointBase::setReferenceCoordSystem(
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem,
		bool theKeepGlobalCoordsSame)
	{
		// No need to update coord data
		bool checkCoordSystem{};
		if (!theKeepGlobalCoordsSame) checkCoordSystem = true;
		else if (theReferenceCoordSystem->isGlobal()) {
			if (c_referenceCoordSystem->isGlobal())
			{
				checkCoordSystem = true;
			}
		}
		else if (theReferenceCoordSystem->equals(*c_referenceCoordSystem))
		{
			checkCoordSystem = true;
		}
		if (checkCoordSystem) {
			ReferenceObject::setReferenceCoordSystemBase(theReferenceCoordSystem);
			return;
		}

		// The localcoords need update
		if (theReferenceCoordSystem->isGlobal())
		{
			setLocalCoords(getGlobalCoords());
		}
		else
		{
			Point3D point{ getGlobalCoords() };
			setLocalCoords(theReferenceCoordSystem->measurePointCoords(point));
		}

		// Set the reference CS
		ReferenceObject::setReferenceCoordSystemBase(theReferenceCoordSystem);
	}

	/// <summary>
	/// Setter - local coord - X
	/// </summary>
	void PointBase::setLocalCoordX(const double& theLocalCoordX)
	{
		c_localCoords[0] = theLocalCoordX;
	}

	/// <summary>
	/// Setter - local coord - Y
	/// </summary>
	void PointBase::setLocalCoordY(const double& theLocalCoordY)
	{
		c_localCoords[1] = theLocalCoordY;
	}

	/// <summary>
	/// Setter - local coord - Z
	/// </summary>
	void PointBase::setLocalCoordZ(const double& theLocalCoordZ)
	{
		c_localCoords[2] = theLocalCoordZ;
	}

	/// <summary>
	/// Setter - local coords
	/// </summary>
	void PointBase::setLocalCoords(const std::array<double, 3>& theLocalCoords)
	{
		c_localCoords = theLocalCoords;
	}

	/// <summary>
	/// Setter - local coords
	/// </summary>
	void PointBase::setLocalCoords(const std::vector<double, std::allocator<double>>& theLocalCoords)
	{
		std::copy(theLocalCoords.begin(), theLocalCoords.end(), c_localCoords.begin());
	}

	/// <summary>
	/// Returns if the point coincides with the input point
	/// </summary>
	/// <exception> NullptrException </exception>
	bool PointBase::coincides(ARGCOPY(PointBase) thePoint) const
	{
		// Get the item in my reference CS
		return GeometryMath::equals(getGlobalCoords(), thePoint.getGlobalCoords(), getToleranceGeneral());
	}

	/// <summary>
	/// Calculates distance to point
	/// </summary>
	/// <exception> NullptrException </exception>
	double PointBase::calculateDistance(ARGCOPY(PointBase) thePoint) const
	{
		auto globalCoords0{ getGlobalCoords() };
		auto globalCoords1{ thePoint.getGlobalCoords() };
		double outDistance{};
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++)
		{
			outDistance += std::pow(globalCoords0[iCoord] - globalCoords1[iCoord], 2.);
		}
		return std::pow(outDistance, 0.5);
	}

	/// <summary>
	/// Creates mid-point
	/// </summary>
	auto PointBase::createMidPoint(ARGCOPY(PointBase) thePoint) const -> std::shared_ptr<PointBase>
	{
		return createInterpolationPoint(thePoint, 0.5);
	}

	/// <summary>
	/// Creates interpolation point
	/// theFactor: How much of the distance between the two values is to be travelled
	///	Betwen 0. and 1.
	/// </summary>
	/// <exception> NullptrException </exception>
	auto PointBase::createInterpolationPoint(
		ARGCOPY(PointBase) thePoint,
		const double& theFactor) const
		-> std::shared_ptr<PointBase>
	{
		// Get the item in my reference CS
		auto point = getPointWithMyCoordSystem(thePoint);
		auto interpolationCoords = PointBase::interpolateCoords(c_localCoords, point->getLocalCoords(), theFactor);
		if (is2D() && point->is2D())
		{
			return std::make_shared<Point2D>(
				c_referenceCoordSystem,
				interpolationCoords);
		}
		return std::make_shared<Point3D>(
			c_referenceCoordSystem,
			interpolationCoords);
	}

	/// <summary>
	/// Creates interpolation point
	/// theFactor: How much of the distance between the two values is to be travelled
	///	Betwen 0. and 1.
	/// Reference CS is the global CS
	/// </summary>
	std::array<double, 3> PointBase::interpolateCoords(
		const std::array<double, 3>& theCoords0,
		const std::array<double, 3>& theCoords1,
		const double& theFactor)
	{
		std::array<double, 3> outCoords = { {} };
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			outCoords[iCoord] = GeometryMath::interpolate(theCoords0[iCoord], theCoords1[iCoord], theFactor);
		}
		return outCoords;
	}

	/// <summary>
	/// Creates point at the origin of the global CS
	/// </summary>
	auto PointBase::createPointAtOrigin(int theDimensionCount) -> std::shared_ptr<PointBase>
	{
		if (theDimensionCount == DIMENSIONS::D2)
		{
			return std::make_shared<Point2D>(std::array<double, 3>{{}});
		}
		return std::make_shared<Point3D>(std::array<double, 3>{{}});
	}

	/// <summary>
	/// Creates point at the origin of the input CS
	/// </summary>
	auto PointBase::createPointAtOrigin(
		const int theDimensionCount,
		const std::shared_ptr<CoordSystem>& theReferenceCoordSystem)
		-> std::shared_ptr<PointBase>
	{
		if (theDimensionCount == DIMENSIONS::D2)
		{
			return std::make_shared<Point2D>(
				theReferenceCoordSystem,
				std::array<double, 3>{{}});
		}
		return std::make_shared<Point3D>(
			theReferenceCoordSystem,
			std::array<double, 3>{{}});
	}

	/// <summary>
	/// See GeometryObject::Clone(const T& arg) for descriptions.
	/// </summary>
	template<typename T>
	T PointBase::Clone(const T& arg)
	{
		T point = GeometryObject::Clone(arg);
		point.setMembers(
			arg.getDimensionCount(),
			arg.getReferenceCoordSystem(),
			arg.getLocalCoords());
		return point;
	}
}
