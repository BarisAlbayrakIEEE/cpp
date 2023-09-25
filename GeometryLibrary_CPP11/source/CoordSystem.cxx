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
	/// Protected ctor which acts like a setter
	/// Preconditions are not inspected (e.g. is X perpandicular to Y and Z?)
	/// An internal ctor
	/// </summary>
	CoordSystem::CoordSystem(
		const std::array<double, 3>& theOriginCoords,
		const std::array<double, 3>& theAxisComponentsX,
		const std::array<double, 3>& theAxisComponentsY,
		const std::array<double, 3>& theAxisComponentsZ)
	:
	GeometryObject(TOLERANCE_GENERAL, TOLERANCE_SENSITIVE),
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
	/// Set the members.
	/// 
	/// Follows RAII idiom
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
	/// Set the members.
	/// 
	/// Follows RAII idiom
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
	/// This operator inspects direct equality which requires direct equality of all members.
	/// The CoordSystem class does not have geometrical equality/unequality methods separately
	/// (e.g. += operator of PointBase class)
	/// as the members of CoordSystem class are of fundamental type (array)
	/// </summary>
	bool CoordSystem::operator==(const CoordSystem& rhs) const
	{
		return equals(rhs);
	}

	/// <summary>
	/// See the docstring of == operator
	/// </summary>
	bool CoordSystem::operator!=(const CoordSystem& rhs) const
	{
		return !operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for CoordSystem.
	/// Hence, same as == operator
	/// </summary>
	bool CoordSystem::operator+=(const CoordSystem& rhs) const
	{
		return operator==(rhs);
	}

	/// <summary>
	/// Equality and geometrical equality are the same for CoordSystem.
	/// Hence, same as != operator
	/// </summary>
	bool CoordSystem::operator-=(const CoordSystem& rhs) const
	{
		return !operator+=(rhs);
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
	/// This method inspects direct equality which requires direct equality of all members.
	/// </summary>
	bool CoordSystem::equals(ARGCOPY(CoordSystem) theCoordSystem) const {
		if (this == &theCoordSystem) return true;
		if (c_isGlobal != theCoordSystem.c_isGlobal) return false;
		if (!GeometryMath::equals(c_originCoords, theCoordSystem.getOriginCoords(), getToleranceGeneral())) return false;
		if (!GeometryMath::equals(c_axisComponentsX, theCoordSystem.getAxisComponentsX(), getToleranceGeneral())) return false;
		if (!GeometryMath::equals(c_axisComponentsY, theCoordSystem.getAxisComponentsY(), getToleranceGeneral())) return false;
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
	auto CoordSystem::getOriginCoords() const -> std::array<double, 3>
	{
		return c_originCoords;
	}

	/// <summary>
	/// Getter - The components of the axis - X
	/// </summary>
	auto CoordSystem::getAxisComponentsX() const -> std::array<double, 3>
	{
		return c_axisComponentsX;
	}

	/// <summary>
	/// Getter - The components of the axis - Y
	/// </summary>
	auto CoordSystem::getAxisComponentsY() const -> std::array<double, 3>
	{
		return c_axisComponentsY;
	}

	/// </summary>
	/// Getter - The components of the axis - Z
	/// </summary>
	auto CoordSystem::getAxisComponentsZ() const -> std::array<double, 3>
	{
		return c_axisComponentsZ;
	}

	/// </summary>
	/// Getter - The axis as a vector - X
	/// </summary>
	auto CoordSystem::getAxisAsVectorX() const -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(c_axisComponentsX);
	}

	/// </summary>
	/// Getter - The axis as a vector - Y
	/// </summary>
	auto CoordSystem::getAxisAsVectorY() const -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(c_axisComponentsY);
	}

	/// </summary>
	/// Getter - The axis as a vector - Z
	/// </summary>
	auto CoordSystem::getAxisAsVectorZ() const -> std::shared_ptr<Vector3D>
	{
		return std::make_shared<Vector3D>(c_axisComponentsZ);
	}

	/// </summary>
	/// Getter - All axes - As Vector - In an array
	/// </summary>
	auto CoordSystem::getAxesAsVector() const -> std::vector<std::shared_ptr<Vector3D>>
	{
		std::vector<std::shared_ptr<Vector3D>> axesAsVector;
		axesAsVector.push_back(getAxisAsVectorX());
		axesAsVector.push_back(getAxisAsVectorY());
		axesAsVector.push_back(getAxisAsVectorZ());
		return axesAsVector;
	}

	/// </summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void CoordSystem::setMembers(
		ARGCOPY(Point3D) theOriginPoint,
		ARGCOPY(Point3D) thePointOnAxisX,
		ARGCOPY(Point3D) thePointOnAxisY)
	{
		Vector3D axisVectorX{ theOriginPoint, thePointOnAxisX };
		Vector3D axisVectorY{ theOriginPoint, thePointOnAxisY };
		setMembers(theOriginPoint, axisVectorX, axisVectorY);
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

		auto axisVectorZ = theAxisVectorX.crossProduct(theAxisVectorY);
		auto axisVectorY = axisVectorZ->crossProduct(theAxisVectorX);
		auto axisX{ GeometryMath::factorizeArray(theAxisVectorX.getGlobalComponents(), theAxisVectorX.getMagnitude()) };
		auto axisY{ GeometryMath::factorizeArray(axisVectorY->getGlobalComponents(), theAxisVectorX.getMagnitude()) };
		auto axisZ{ GeometryMath::factorizeArray(axisVectorZ->getGlobalComponents(), theAxisVectorX.getMagnitude()) };

		setMembers(theOriginPoint.getGlobalCoords(), axisX, axisY, axisZ);
	}

	/// </summary>
	/// Setter - members
	/// Protected method used in this class and child classes only.
	/// </summary>
	void CoordSystem::setMembers(
		const std::array<double, 3>& theOriginCoords,
		const std::array<double, 3>& theAxisComponentsX,
		const std::array<double, 3>& theAxisComponentsY,
		const std::array<double, 3>& theAxisComponentsZ)
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
	void CoordSystem::setOriginCoords(const std::array<double, 3>& theOriginCoords) {
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
		if (!GeometryMath::equals(c_axisComponentsX, theCoordSystem.getAxisComponentsX(), getToleranceGeneral())) return false;
		if (!GeometryMath::equals(c_axisComponentsY, theCoordSystem.getAxisComponentsY(), getToleranceGeneral())) return false;
		return true;
	}

	/// <summary>
	/// Determines if the CSs are identical
	/// </summary>
	bool CoordSystem::isIdentical(ARGCOPY(CoordSystem) theCoordSystem) const
	{
		if (c_isGlobal != theCoordSystem.c_isGlobal) return false;
		if (!GeometryMath::equals(c_originCoords, theCoordSystem.getOriginCoords(), getToleranceGeneral())) return false;
		return isParallel(theCoordSystem);
	}

	/// <summary>
	/// Creates a point having this CS as the reference CS and with the input coords
	/// </summary>
	auto CoordSystem::createPoint(const std::array<double, 3>& theCoords) {
		return std::make_shared<Point3D>(std::shared_ptr<CoordSystem>(this), theCoords);
	}

	/// <summary>
	/// Measures coords of the input point wrt this CS.
	/// Returns the local coords of the point if the reference CS of the point is this CS
	/// </summary>
	auto CoordSystem::measurePointCoords(ARGCOPY(PointBase) thePoint) const -> std::array<double, 3>
	{
		// Check if the point reference CS is the same
		if (equals(*thePoint.getReferenceCoordSystem())) return thePoint.getLocalCoords();

		// Determine the inverse of the matrix
		// No need to catch exception as not possible (member inspection performed at the beginning)
		std::array<std::array<double, 3>, 3> basisVectors = {
			std::array<double, 3>{{}},
			std::array<double, 3>{{}},
			std::array<double, 3>{{}} };
		basisVectors[0] = c_axisComponentsX;
		basisVectors[1] = c_axisComponentsY;
		basisVectors[2] = c_axisComponentsZ;
		auto inverseMatrix = GeometryMath::calculateMatrixInverseS33(basisVectors, getToleranceSensitive());

		// Measure the point coords wrt this
		auto originCoords = getOriginCoords();
		auto globalCoords = thePoint.getGlobalCoords();
		std::array<double, 3> outCoords = { {} };
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
	auto CoordSystem::measureVectorComponents(ARGCOPY(VectorBase) theVector) const -> std::array<double, 3>
	{
		// Check if the vector reference CS is the same
		if (equals(*theVector.getReferenceCoordSystem())) return theVector.getLocalComponents();

		// Determine the inverse of the matrix
		// No need to catch exception as not possible (member inspection performed at the beginning)
		std::array<std::array<double, 3>, 3> basisVectors = {
			std::array<double, 3>{{}},
			std::array<double, 3>{{}},
			std::array<double, 3>{{}} };
		basisVectors[0] = c_axisComponentsX;
		basisVectors[1] = c_axisComponentsY;
		basisVectors[2] = c_axisComponentsZ;
		auto inverseMatrix = GeometryMath::calculateMatrixInverseS33(basisVectors, getToleranceSensitive());

		return GeometryMath::multiplyMatrixToVector<double, 3>(inverseMatrix, theVector.getGlobalComponents());
	}

	/// <summary>
	/// Creates a new point by rotating the input point about x-axis
	/// </summary>
	auto CoordSystem::rotatePointAboutAxisX(ARGCOPY(PointBase) thePoint, double theAngle) {
		return rotatePointBase(thePoint, theAngle, 1, 2, 0);
	}

	/// <summary>
	/// Creates a new point by rotating the input point about y-axis
	/// </summary>
	auto CoordSystem::rotatePointAboutAxisY(ARGCOPY(PointBase) thePoint, double theAngle) {
		return rotatePointBase(thePoint, theAngle, 2, 0, 1);
	}

	/// <summary>
	/// Creates a new point by rotating the input point about z-axis
	/// </summary>
	auto CoordSystem::rotatePointAboutAxisZ(ARGCOPY(PointBase) thePoint, double theAngle) {
		return rotatePointBase(thePoint, theAngle, 0, 1, 2);
	}

	/// <summary>
	/// Method:
	///	The rotation will be performed about an axis
	///	The coordinate about the reference axis (theAxis2 input) remains the same for this operation
	///	The other two (theAxis0 and theAxis1) are simply rotated
	/// theAxis0, theAxis1 and theAxis2 inputs are the indicies of the axes: One of [0,3)
	/// </summary>
	auto CoordSystem::rotatePointBase(
		ARGCOPY(PointBase) thePoint,
		const double& theAngle,
		int theAxis0,
		int theAxis1,
		int theAxis2) -> std::shared_ptr<Point3D>
	{
		// Calculate the radius and the current angle assumming that the rotation traces a circle centered at the origin
		auto globalCoords = thePoint.getGlobalCoords();
		double angleCurrent = std::atan(globalCoords[theAxis1] / globalCoords[theAxis0]);
		double radius = std::pow(std::pow(globalCoords[theAxis0], 2.) + std::pow(globalCoords[theAxis1], 2.), 0.5);

		// Calculate the final values for the rotated coordinates
		std::array<double, 3> rotatedCoords{
			radius * std::cos(angleCurrent + theAngle),
			radius * std::sin(angleCurrent + theAngle),
			globalCoords[theAxis2] };
		return std::make_shared<Point3D>(std::shared_ptr<CoordSystem>(this), rotatedCoords);
	}
}
