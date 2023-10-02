/// <summary>
/// Includes some simple math functions.
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#include "GeometryMath.hxx"

namespace GeometryNamespace {
	/// <summary/>
	/// Calculates the determinant of a 3x3 matrix defined by a nested std::array
	/// </summary>
	auto GeometryMath::calculateMatrixDeterminantS33(const std::array<std::array<double, 3>, 3>& theMatrix) -> double
	{
		double outDeterminant = 0.;
		outDeterminant += theMatrix[0][0] * (theMatrix[1][1] * theMatrix[2][2] - theMatrix[1][2] * theMatrix[2][1]);
		outDeterminant -= theMatrix[0][1] * (theMatrix[1][0] * theMatrix[2][2] - theMatrix[1][2] * theMatrix[2][0]);
		outDeterminant += theMatrix[0][2] * (theMatrix[1][0] * theMatrix[2][1] - theMatrix[1][1] * theMatrix[2][0]);
		return outDeterminant;
	}

	/// <summary/>
	/// Calculates the inverse of a 3x3 matrix defined by a nested std::array
	/// </summary>
	auto GeometryMath::calculateMatrixInverseS33(
		const std::array<std::array<double, 3>, 3>& theMatrix)
		->std::array<std::array<double, 3>, 3>
	{
		// Determine determinant
		double determinant = calculateMatrixDeterminantS33(theMatrix);

		// Check if inversible matrix
		if (zero_g(determinant)) throw ZeroDeterminantException();

		// The inverse matrix
		std::array<std::array<double, 3>, 3> outInverseMatrix;
		outInverseMatrix[0][0] = (theMatrix[1][1] * theMatrix[2][2] - theMatrix[1][2] * theMatrix[2][1]) / determinant;
		outInverseMatrix[0][1] = (theMatrix[0][2] * theMatrix[2][1] - theMatrix[0][1] * theMatrix[2][2]) / determinant;
		outInverseMatrix[0][2] = (theMatrix[0][1] * theMatrix[1][2] - theMatrix[0][2] * theMatrix[1][1]) / determinant;

		outInverseMatrix[1][0] = (theMatrix[1][2] * theMatrix[2][0] - theMatrix[1][0] * theMatrix[2][2]) / determinant;
		outInverseMatrix[1][1] = (theMatrix[0][0] * theMatrix[2][2] - theMatrix[0][2] * theMatrix[2][0]) / determinant;
		outInverseMatrix[1][2] = (theMatrix[0][2] * theMatrix[1][0] - theMatrix[0][0] * theMatrix[1][2]) / determinant;

		outInverseMatrix[2][0] = (theMatrix[1][0] * theMatrix[2][1] - theMatrix[1][1] * theMatrix[2][0]) / determinant;
		outInverseMatrix[2][1] = (theMatrix[0][1] * theMatrix[2][0] - theMatrix[0][0] * theMatrix[2][1]) / determinant;
		outInverseMatrix[2][2] = (theMatrix[0][0] * theMatrix[1][1] - theMatrix[0][1] * theMatrix[1][0]) / determinant;

		return outInverseMatrix;
	}
}
