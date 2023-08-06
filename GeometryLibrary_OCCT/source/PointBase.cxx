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
	IMPLEMENT_STANDARD_RTTIEXT(PointBase, ReferenceObject)

	/// <summary>
	/// The default constructor
	/// Reference CS is the global CS
	/// Point coords are all default (i.e. point is at the origin of the reference CS)
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(const int theDimensionCount)
		:
		ReferenceObject(theDimensionCount),
		c_localCoords{ {} } { }

	/// <summary>
	/// Ctor
	/// Point by coords
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem) throw (NullptrException)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem)
	{
		c_localCoords = arrayS3{};
	}

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theLocalCoords)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem), c_localCoords{ theLocalCoords } { }

	/// <summary>
	/// The main constructor
	/// CAUTION: Member initialization is not performed to follow RAII
	/// </summary>
	PointBase::PointBase(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theLocalCoords)
		: ReferenceObject(theDimensionCount, theReferenceCoordSystem)
	{
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	PointBase::PointBase(const PointBase& rhs)
		:
		ReferenceObject(DIMENSIONS::D2)
	{
		PointBase::copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	PointBase& PointBase::operator=(const PointBase& rhs) {
		if (&rhs == this) return *this;

		PointBase::copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	PointBase::PointBase(PointBase&& rhs) noexcept
		:
		ReferenceObject(DIMENSIONS::D2)
	{
		PointBase::copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	PointBase& PointBase::operator=(PointBase&& rhs) noexcept
	{
		if (&rhs == this) return *this;

		PointBase::copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The += operator inspects geometrical equality.
	/// </summary>
	bool PointBase::operator==(const PointBase& rhs) {
		if (&rhs == this) return true;
		if (ReferenceObject::operator!=(rhs)) return false;
		return GeometryMath::equals(c_localCoords, rhs.getLocalCoords(), c_toleranceGeneral);
	}

	/// <summary>
	/// This operator inspects direct unequality which requires direct unequality of any member.
	/// The -= operator inspects geometrical unequality.
	/// </summary>
	bool PointBase::operator!=(const PointBase& rhs) {
		return !operator==(rhs);
	}

	/// <summary>
	/// This method inspects final geometrical equality which is actually the coincicience.
	/// Point coordinates and vector components are inspected wrt the global CS.
	/// Additionally, inclusion is used rather than the equivalence for the passing points.
	/// </summary>
	bool PointBase::operator+=(const PointBase& rhs) {
		if (&rhs == this) return true;
		return GeometryMath::equals(getGlobalCoords(), rhs.getGlobalCoords(), c_toleranceGeneral);
	}

	/// <summary>
	/// This method inspects final geometrical unequality.
	/// See += operator docstring for the details.
	/// </summary>
	bool PointBase::operator-=(const PointBase& rhs) {
		return !operator+=(rhs);
	}

	/// <summary>
	/// Use OCCT approach
	/// </summary>
	PointBase::~PointBase() {
		Destroy();
	}

	/// <summary>
	/// Use Nullify method of the OCCT Standard_Handle for the object destruction
	/// </summary>
	void PointBase::Destroy() {
		c_localCoords = arrayS3{ 0., 0., 0. };
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void PointBase::copyBase(const PointBase& rhs)
	{
		ReferenceObject::copyBase(rhs);
		std::copy(std::begin(rhs.getLocalCoords()), std::end(rhs.getLocalCoords()), std::begin(c_localCoords));
	}

	/// <summary>
	/// Base method for both the direct equality and the geometrical equality
	/// </summary>
	bool PointBase::equalsBase(ARGCOPY(PointBase) thePoint) const {
		if (thePoint.IsNull()) return false;
		return true;
	}

	/// <summary>
	/// See == operator docstring
	/// </summary>
	bool PointBase::equals(ARGCOPY(PointBase) thePoint) const {
		if (this == thePoint.get()) return true;
		if (!ReferenceObject::equals(thePoint)) return false;
		if (!equalsBase(thePoint)) return false;
		return GeometryMath::equals(c_localCoords, thePoint->getLocalCoords(), c_toleranceGeneral);
	}

	/// <summary>
	/// See += operator docstring
	/// </summary>
	bool PointBase::equalsGeometrically(ARGCOPY(PointBase) thePoint) const {
		if (this == thePoint.get()) return true;
		if (!equalsBase(thePoint)) return false;
		Handle(PointBase) point = getPointWithMyCoordSystem(thePoint);
		return GeometryMath::equals(c_localCoords, point->getLocalCoords(), c_toleranceGeneral);
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
	arrayS3 PointBase::getLocalCoords() const {
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
	double PointBase::getGlobalCoordZ() const {
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
	arrayS3 PointBase::getGlobalCoords() const {
		if (c_referenceCoordSystem.IsNull()) return c_localCoords;
		if (c_referenceCoordSystem->isGlobal()) return c_localCoords;

		std::vector<Handle(Vector3D)> axesVectors{ c_referenceCoordSystem->getAxesAsVector() }; // By definition, the reference CS is the global CS
		arrayS3 outGlobalCoords{ c_referenceCoordSystem->getOriginCoords() };
		for (int iAxis = 0; iAxis < DIMENSIONS::D3; iAxis++) {
			arrayS3 vectorComponents{ axesVectors[iAxis]->getLocalComponents() };
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
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const arrayS3& theLocalCoords)
	{
		ReferenceObject::setMembers(theDimensionCount, theReferenceCoordSystem);
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void PointBase::setMembers(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
		const vectorInput1D& theLocalCoords)
	{
		ReferenceObject::setMembers(theDimensionCount, theReferenceCoordSystem);
		setLocalCoords(theLocalCoords);
	}

	/// <summary>
	/// Setter - reference CS
	/// Public method for the users.
	/// theKeepGlobalCoordsSarne:
	///		true:
	///			Keep the coords wrt the GLOBAL CS same and update the coords wrt the REFERENCE CS.
	///		false:
	///			Keep the coords wrt the REFERENCE CS same and update the coords wrt the GLOBAL CS.
	/// </summary>
	void PointBase::setReferenceCoordSystem(
		ARGCOPY(CoordSystem) theReferenceCoordSystem,
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
		else if (theReferenceCoordSystem->equals(c_referenceCoordSystem))
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
			Handle(Point3D) point{ new Point3D(GlobalCoordSystem::getGlobalCoordSystem(), getGlobalCoords()) };
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
	void PointBase::setLocalCoords(const arrayS3& theLocalCoords)
	{
		c_localCoords = theLocalCoords;
	}

	/// <summary>
	/// Setter - local coords
	/// </summary>
	void PointBase::setLocalCoords(const vectorInput1D& theLocalCoords)
	{
		inspectVectorInput(theLocalCoords);
		c_localCoords = GeometryMath::convertVectorToArray1DS3(theLocalCoords);
	}

	/// <summary>
	/// Returns if the point coincides with the input point
	/// </summary>
	/// <exception> NullptrException </exception>
	bool PointBase::coincides(ARGCOPY(PointBase) thePoint) const throw (NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Get the item in my reference CS
		return GeometryMath::equals(getGlobalCoords(), thePoint->getGlobalCoords(), c_toleranceGeneral);
	}

	/// <summary>
	/// Calculates distance to point
	/// </summary>
	/// <exception> NullptrException </exception>
	double PointBase::calculateDistance(ARGCOPY(PointBase) thePoint) const throw (NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		arrayS3 globalCoords0{ getGlobalCoords() };
		arrayS3 globalCoords1{ thePoint->getGlobalCoords() };
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
	OUTVAL(PointBase) PointBase::createMidPoint(ARGCOPY(PointBase) thePoint) const {
		return createInterpolationPoint(thePoint, 0.5);
	}

	/// <summary>
	/// Creates interpolation point
	/// theFactor: How much of the distance between the two values is to be travelled
	///	Betwen 0. and 1.
	/// </summary>
	/// <exception> NullptrException </exception>
	OUTVAL(PointBase) PointBase::createInterpolationPoint(
		ARGCOPY(PointBase) thePoint,
		const double& theFactor) const throw (NullptrException)
	{
		if (thePoint.IsNull())
		{
			throw NullptrException();
		}

		// Get the item in my reference CS
		Handle(PointBase) point = getPointWithMyCoordSystem(thePoint);

		arrayS3 interpolationCoords = PointBase::interpolateCoords(c_localCoords, point->getLocalCoords(), theFactor);
		if (is2D() && point->is2D())
		{
			return Handle(Point2D)(new Point2D(c_referenceCoordSystem, interpolationCoords));
		}
		else
		{
			return Handle(Point3D)(new Point3D(c_referenceCoordSystem, interpolationCoords));
		}
	}

	/// <summary>
	/// Creates interpolation point
	/// theFactor: How much of the distance between the two values is to be travelled
	///	Betwen 0. and 1.
	/// Reference CS is the global CS
	/// </summary>
	arrayS3 PointBase::interpolateCoords(const arrayS3& theCoords0, const arrayS3& theCoords1, const double& theFactor)
	{
		arrayS3 outCoords;
		for (int iCoord = 0; iCoord < DIMENSIONS::D3; iCoord++) {
			outCoords[iCoord] = GeometryMath::interpolate(theCoords0[iCoord], theCoords1[iCoord], theFactor);
		}
		return outCoords;
	}

	/// <summary>
	/// Creates point at the origin of the global CS
	/// </summary>
	OUTVAL(PointBase) PointBase::createPointAtOrigin(int theDimensionCount)
	{
		if (theDimensionCount == DIMENSIONS::D2)
		{
			return Handle(Point2D)(new Point2D(GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{}));
		}
		else
		{
			return Handle(Point3D)(new Point3D(GlobalCoordSystem::getGlobalCoordSystem(), arrayS3{}));
		}
	}

	/// <summary>
	/// Creates point at the origin of the input CS
	/// </summary>
	OUTVAL(PointBase) PointBase::createPointAtOrigin(
		const int theDimensionCount,
		ARGCOPY(CoordSystem) theReferenceCoordSystem)
	{
		if (theDimensionCount == DIMENSIONS::D2)
		{
			return Handle(Point2D)(new Point2D(theReferenceCoordSystem, arrayS3{}));
		}
		else
		{
			return Handle(Point3D)(new Point3D(theReferenceCoordSystem, arrayS3{}));
		}
	}
}
