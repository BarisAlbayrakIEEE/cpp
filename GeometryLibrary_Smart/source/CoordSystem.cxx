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
	/// Protected ctor to be used in the project which acts like a setter
	/// </summary>
	CoordSystem::CoordSystem(
		const arrayS3& theOriginCoords,
		const arrayS3& theAxisComponentsX,
		const arrayS3& theAxisComponentsY,
		const arrayS3& theAxisComponentsZ)
	:
	GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
	c_isGlobal{ false },
	c_originCoords{ theOriginCoords },
	c_axisComponentsX{ theAxisComponentsX },
	c_axisComponentsY{ theAxisComponentsY },
	c_axisComponentsZ{ theAxisComponentsZ } { }

	/// <summary>
	/// Ctor with three points:
	/// Create two vectors simulating axes x and y.
	/// However, y-axis will be updated later as the two vectors may not be perpandicular.
	/// Axis z is obtained using the cross product of the vectors: Z = X x Y.
	/// Recalculate y-axis using the cross product of the vectors: Y = Z x X
	/// Set the membbers.
	/// </summary>
	CoordSystem::CoordSystem(
		ARGCOPY(Point3D) theOriginPoint,
		ARGCOPY(Point3D) thePointOnAxisX,
		ARGCOPY(Point3D) thePointOnAxisY)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(theOriginPoint, thePointOnAxisX, thePointOnAxisY);
	}

	/// <summary>
	/// Ctor with a point andd two vectors:
	/// y-axis will be updated later as the two vectors may not be perpandicular.
	/// Axis z is obtained using the cross product of the vectors: Z = X x Y.
	/// Recalculate y-axis using the cross product of the vectors: Y = Z x X
	/// Set the membbers.
	/// </summary>
	CoordSystem::CoordSystem(
		ARGCOPY(Point3D) theOriginPoint,
		ARGCOPY(Vector3D) theAxisVectorX,
		ARGCOPY(Vector3D) theAxisVectorY)
		: GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE)
	{
		setMembers(theOriginPoint, theAxisVectorX, theAxisVectorY);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	CoordSystem::CoordSystem(const CoordSystem& rhs)
	{
		copyBase(rhs);
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	CoordSystem& CoordSystem::operator=(const CoordSystem& rhs) {
		if (&rhs == this) return *this;

		copyBase(rhs);
		return *this;
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	CoordSystem::CoordSystem(CoordSystem&& rhs) noexcept
	{
		copyBase(rhs);
		rhs.Destroy();
	}

	/// <summary>
	/// Keep in mind that the copy and move ctors and assignment operators are deleted from
	/// GeometryObject, ReferenceObject and VectorBase
	/// in order to prevent slicing
	/// </summary>
	CoordSystem& CoordSystem::operator=(CoordSystem&& rhs) noexcept {
		if (&rhs == this) return *this;

		copyBase(rhs);
		rhs.Destroy();
		return *this;
	}

	/// <summary>
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The CoordSystem class does not have geometrical equality/unequality methods separately
	/// (e.g. += operator of PointBase class)
	/// as the members of CoordSystem class are of fundamental type (array)
	/// </summary>
	bool CoordSystem::operator==(const CoordSystem& rhs) {
		if (&rhs == this) return true;
		if (!GeometryMath::equals(c_originCoords, rhs.getOriginCoords(), c_toleranceGeneral)) return false;
		if (!GeometryMath::equals(c_axisComponentsX, rhs.getAxisComponentsX(), c_toleranceGeneral)) return false;
		if (!GeometryMath::equals(c_axisComponentsY, rhs.getAxisComponentsY(), c_toleranceGeneral)) return false;
		if (!GeometryMath::equals(c_axisComponentsZ, rhs.getAxisComponentsZ(), c_toleranceGeneral)) return false;
		return true;
	}

	/// <summary>
	/// See the docstring of == operator
	/// </summary>
	bool CoordSystem::operator!=(const CoordSystem& rhs) {
		return !operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for CoordSystem.
	/// Hence, same as == operator
	/// </summary>
	bool CoordSystem::operator+=(const CoordSystem& rhs)
	{
		return operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for CoordSystem.
	/// Hence, same as != operator
	/// </summary>
	bool CoordSystem::operator-=(const CoordSystem& rhs)
	{
		return !operator+=(rhs);
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	CoordSystem::~CoordSystem() {
		Destroy();
	}

	/// <summary>
	/// Used in the copy/move ctor and operators
	/// </summary>
	void CoordSystem::copyBase(const CoordSystem& rhs)
	{
		GeometryObject::copyBase(rhs);
		c_isGlobal = rhs.c_isGlobal;
		std::copy(std::begin(rhs.c_originCoords), std::end(rhs.c_originCoords), std::begin(c_originCoords));
		std::copy(std::begin(rhs.c_axisComponentsX), std::end(rhs.c_axisComponentsX), std::begin(c_axisComponentsX));
		std::copy(std::begin(rhs.c_axisComponentsY), std::end(rhs.c_axisComponentsY), std::begin(c_axisComponentsY));
		std::copy(std::begin(rhs.c_axisComponentsZ), std::end(rhs.c_axisComponentsZ), std::begin(c_axisComponentsZ));
	}

	/// <summary>
	/// Pure virtual function from GeometryAnstractObject
	/// </summary>
	bool CoordSystem::is2D() const
	{
		return false;
	}

	/// <summary>
	/// Pure virtual function from GeometryAnstractObject
	/// </summary>
	bool CoordSystem::is3D() const
	{
		return true;
	}

	/// <summary>
	/// Actually, is the defaault dtor which is not a good approach to explicitly write the default dtor
	/// However, kept explicitly in the code due to the class hierarchy and slicing issue.
	/// </summary>
	void CoordSystem::Destroy() {
		c_originCoords = arrayS3{};
		c_axisComponentsX = arrayS3{};
		c_axisComponentsY = arrayS3{};
		c_axisComponentsZ = arrayS3{};
	}

	/// <summary>
	/// This method inspects direct equality which requires direct equality of all members.
	/// </summary>
	bool CoordSystem::equals(ARGCOPY(CoordSystem) theCoordSystem) const {
		if (this == &theCoordSystem) return true;
		if (c_isGlobal != theCoordSystem.c_isGlobal) return false;
		if (!GeometryMath::equals(c_originCoords, theCoordSystem.getOriginCoords(), c_toleranceGeneral)) return false;
		if (!GeometryMath::equals(c_axisComponentsX, theCoordSystem.getAxisComponentsX(), c_toleranceGeneral)) return false;
		if (!GeometryMath::equals(c_axisComponentsY, theCoordSystem.getAxisComponentsY(), c_toleranceGeneral)) return false;
		if (!GeometryMath::equals(c_axisComponentsZ, theCoordSystem.getAxisComponentsZ(), c_toleranceGeneral)) return false;
		return true;
	}

	/// <summary>
	/// Equality and geometrical equality are the same for CoordSystem.
	/// Hence, same as equals(theCoordSystem)
	/// </summary>
	bool CoordSystem::equalsGeometrically(ARGCOPY(CoordSystem) theCoordSystem) const
	{
		return equals(theCoordSystem);
	}

	/// <summary>
	/// Getter - Origin coords (global coords)
	/// </summary>
	arrayS3 CoordSystem::getOriginCoords() const {
		return c_originCoords;
	}

	/// <summary>
	/// Getter - The components of the axis - X
	/// </summary>
	arrayS3 CoordSystem::getAxisComponentsX() const {
		return c_axisComponentsX;
	}

	/// <summary>
	/// Getter - The components of the axis - Y
	/// </summary>
	arrayS3 CoordSystem::getAxisComponentsY() const {
		return c_axisComponentsY;
	}

	/// </summary>
	/// Getter - The components of the axis - Z
	/// </summary>
	arrayS3 CoordSystem::getAxisComponentsZ() const {
		return c_axisComponentsZ;
	}

	/// </summary>
	/// Getter - The axis as a vector - X
	/// </summary>
	Vector3D CoordSystem::getAxisAsVectorX() const {
		return Vector3D(c_axisComponentsX);
	}

	/// </summary>
	/// Getter - The axis as a vector - Y
	/// </summary>
	Vector3D CoordSystem::getAxisAsVectorY() const {
		return Vector3D(c_axisComponentsY);
	}

	/// </summary>
	/// Getter - The axis as a vector - Z
	/// </summary>
	Vector3D CoordSystem::getAxisAsVectorZ() const {
		return Vector3D(c_axisComponentsZ);
	}

	/// </summary>
	/// Getter - All axes - As Vector - In an array
	/// </summary>
	std::vector<Vector3D> CoordSystem::getAxesAsVector() const {
		std::vector<Vector3D> axesAsVector;
		axesAsVector.push_back(getAxisAsVectorX());
		axesAsVector.push_back(getAxisAsVectorY());
		axesAsVector.push_back(getAxisAsVectorZ());
		return axesAsVector;
	}

	/// </summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </sumrnary>
	void CoordSystem::setMembers(
		ARGCOPY(Point3D) theOriginPoint,
		ARGCOPY(Point3D) thePointOnAxisX,
		ARGCOPY(Point3D) thePointOnAxisY)
	{
		try {
			Vector3D axisVectorX { theOriginPoint, thePointOnAxisX };
			Vector3D axisVectorY { theOriginPoint, thePointOnAxisY };
			setMembers(theOriginPoint, axisVectorX, axisVectorY);
		}
		catch (GeometryException) {
			throw;
		}
	}

	/// <summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void CoordSystem::setMembers(
		ARGCOPY(Point3D) theOriginPoint,
		ARGCOPY(Vector3D) theAxisVectorX,
		ARGCOPY(Vector3D) theAxisVectorY)
	{
		if (theAxisVectorX.isParallel(theAxisVectorY)) throw ParallelAxisException();

		Vector3D axisVectorZ = theAxisVectorX.crossProduct(theAxisVectorY);
		Vector3D axisVectorY = axisVectorZ.crossProduct(theAxisVectorX);
		arrayS3 axisX = GeometryMath::factorizeArray(
			theAxisVectorX.getGlobalComponents(),
			1. / theAxisVectorX.getMagnitude());
		arrayS3 axisY = GeometryMath::factorizeArray(
			axisVectorY.getGlobalComponents(),
			1. / axisVectorY.getMagnitude());
		arrayS3 axisZ = GeometryMath::factorizeArray(
			axisVectorZ.getGlobalComponents(),
			1. / axisVectorZ.getMagnitude());

		setMembers(theOriginPoint.getGlobalCoords(), axisX, axisY, axisZ);
	}

	/// </summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void CoordSystem::setMembers(
		const arrayS3& theOriginCoords,
		const arrayS3& theAxisComponentsX,
		const arrayS3& theAxisComponentsY,
		const arrayS3& theAxisComponentsZ)
	{
		c_originCoords = theOriginCoords;
		c_axisComponentsX = theAxisComponentsX;
		c_axisComponentsY = theAxisComponentsY;
		c_axisComponentsZ = theAxisComponentsZ;
		c_isGlobal = false;
	}

	/// </summary>
	/// Setter - Origin coords (global coords)
	/// </summary>
	void CoordSystem::setOriginCoords(const arrayS3& theOriginCoords) {
		c_originCoords = theOriginCoords;
	}

	/// </summary>
	/// Returns if the CS is global
	/// </summary>
	bool CoordSystem::isGlobal() const {
		return c_isGlobal;
	}

	/// </summary>
	/// Determines if the CSs are parallel
	/// </summary>
	bool CoordSystem::isParallel(ARGCOPY(CoordSystem) theCoordSystem) const
	{
		if (!(getAxisAsVectorX().isInTheSameDirection(theCoordSystem.getAxisAsVectorX())))
		{
			return false;
		}
		if (!(getAxisAsVectorY().isInTheSameDirection(theCoordSystem.getAxisAsVectorY())))
		{
			return false;
		}
		if (!(getAxisAsVectorZ().isInTheSameDirection(theCoordSystem.getAxisAsVectorZ())))
		{
			return false;
		}
		return true;
	}

	/// <summary>
	/// Determines if the CSs are identical
	/// </summary>
	bool CoordSystem::isIdentical(ARGCOPY(CoordSystem) theCoordSystem) const
	{
		if (c_isGlobal != theCoordSystem.c_isGlobal) return false;
		if (!GeometryMath::equals(c_originCoords, theCoordSystem.getOriginCoords(), c_toleranceGeneral)) return false;
		return isParallel(theCoordSystem);
	}

	/// <summary>
	/// Creates a point having this CS as the reference CS and with the input coords
	/// </summary>
	PointBase CoordSystem::createPoint(int theDimensionCount, const arrayS3& theCoords) {
		if (theDimensionCount == DIMENSIONS::D2)
		{
			return Point2D(*this, theCoords);
		}
		return Point3D(*this, theCoords);
	}

	/// <summary>
	/// Measures coords of the input point wrt this CS.
	/// Returns the local coords of the point if the reference CS of the point is this CS
	/// </summary>
	arrayS3 CoordSystem::measurePointCoords(ARGCOPY(PointBase) thePoint) const
	{
		// Check if the point reference CS is the same
		if (equals(thePoint.getReferenceCoordSystem())) return thePoint.getLocalCoords();

		// Determine the inverse of the matrix
		// No need to catch exception as not possible (member inspection performed at the beginning)
		arrayS33 basisVectors;
		basisVectors[0] = c_axisComponentsX;
		basisVectors[1] = c_axisComponentsY;
		basisVectors[2] = c_axisComponentsZ;
		arrayS33 inverseMatrix = GeometryMath::calculateMatrixInverseS33(basisVectors, c_toleranceSensitive);

		// Measure the point coords wrt this
		arrayS3 originCoords = getOriginCoords();
		arrayS3 globalCoords = thePoint.getGlobalCoords();
		arrayS3 outCoords;
		for (int iCoord0 = 0; iCoord0 < DIMENSIONS::D3; iCoord0++) {
			outCoords[iCoord0] = 0.;
			for (int iCoord1 = 0; iCoord1 < DIMENSIONS::D3; iCoord1++)
				outCoords[iCoord0] += inverseMatrix[iCoord0][iCoord1] * (globalCoords[iCoord1] - originCoords[iCoord1]);
		}

		return outCoords;
	}

	/// <summary>
	/// Measures components of the input vector wrt this CS.
	/// Returns the local components of the vector if the reference CS of the vector is this CS
	/// </summary>
	arrayS3 CoordSystem::measureVectorComponents(ARGCOPY(VectorBase) theVector) const
	{
		// Check if the vector reference CS is the same
		if (equals(theVector.getReferenceCoordSystem())) return theVector.getLocalComponents();

		// Determine the inverse of the matrix
		// No need to catch exception as not possible (member inspection performed at th beginning)
		arrayS33 basisVectors;
		basisVectors[0] = c_axisComponentsX;
		basisVectors[1] = c_axisComponentsY;
		basisVectors[2] = c_axisComponentsZ;
		arrayS33 inverseMatrix = GeometryMath::calculateMatrixInverseS33(basisVectors, c_toleranceSensitive);

		return GeometryMath::multiplyMatrixToVectorS33S3(basisVectors, theVector.getGlobalComponents());
	}

	/// <summary>
	/// Creates a new point by rotating the input point about x-axis
	/// </summary>
	Point3D CoordSystem::rotatePointAboutAxisX(ARGCOPY(PointBase) thePoint, double theAngle) {
		return rotatePointBase(thePoint, theAngle, 1, 2, 0);
	}

	/// <summary>
	/// Creates a new point by rotating the input point about y-axis
	/// </summary>
	Point3D CoordSystem::rotatePointAboutAxisY(ARGCOPY(PointBase) thePoint, double theAngle) {
		return rotatePointBase(thePoint, theAngle, 2, 0, 1);
	}

	/// <summary>
	/// Creates a new point by rotating the input point about z-axis
	/// </summary>
	Point3D CoordSystem::rotatePointAboutAxisZ(ARGCOPY(PointBase) thePoint, double theAngle) {
		return rotatePointBase(thePoint, theAngle, 0, 1, 2);
	}

	/// <summary>
	/// Method:
	///	The rotation will be performed about an axis
	///	The coordinate about the reference axis (theAxis2 input) remains the same for this operation
	///	The other two (theAxis0 and theAxis1) are simply rotated
	/// theAxis0, theAxis1 and theAxis2 inputs are the indicies of the axes: One of [0,3)
	/// </sumrnary>
	Point3D CoordSystem::rotatePointBase(
		ARGCOPY(PointBase) thePoint,
		const double& theAngle,
		int theAxis0,
		int theAxis1,
		int theAxis2)
	{
		// Calculate the radius and the current angle assumming that the rotation traces a circle centered at the origin
		arrayS3 globalCoords = thePoint.getGlobalCoords();
		double angleCurrent = std::atan(globalCoords[theAxis1] / globalCoords[theAxis0]);
		double radius = std::pow(std::pow(globalCoords[theAxis0], 2.) + std::pow(globalCoords[theAxis1], 2.), 0.5);

		// Calculate the final values for the rotated coordinates
		arrayS3 rotatedCoords{
			radius * std::cos(angleCurrent + theAngle),
			radius * std::sin(angleCurrent + theAngle),
			globalCoords[theAxis2] };
		Point3D point { *this, rotatedCoords };
		return point;
	}
}
