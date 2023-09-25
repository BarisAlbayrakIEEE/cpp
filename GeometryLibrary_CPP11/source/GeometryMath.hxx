/// <summary>
/// Includes some simple math functions.
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#ifndef _GeometryMath_HeaderFile
#define _GeometryMath_HeaderFile

#ifndef _Macros_HeaderFile
#include "Macros.h"
#endif
#ifndef _GeometryException_HeaderFile
#include "GeometryException.hxx"
#endif
#ifndef _GeometryParameters_HeaderFile
#include "GeometryParameters.hxx"
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <array>
#include <vector>
#include <algorithm>

namespace GeometryNamespace {
	class GeometryMath
	{
	public:
		/// <summary>
		///	Interpolates between the two input values
		/// theFactor: Betwen 0. and 1.: [0., 1.]
		/// </summary>
		static double interpolate(const double& theValue0, const double& theValue1, const double& theFactor)
		{
			if (theFactor <= 0.) return theValue0;
			if (theFactor >= 1.) return theValue1;
			return theValue0 + (theValue1 - theValue0) * theFactor;
		}

		/// <summary/>
		/// Inspects equality (using a tolerance value) between two inputs
		/// </summary>
		static bool equals(const double& theValue0, const double& theValue1, const double& theTolerance)
		{
			return std::fabs(theValue0 - theValue1) < theTolerance;
		}

		/// <summary/>
		/// Inspects element-wise equality (using a tolerance value) between two std::arrays
		/// </summary>
		template <typename T, size_t M>
		static bool equals(
			const std::array<T, M>& theValues0,
			const std::array<T, M>& theValues1,
			const double& theTolerance) {
			for (int i = 0; i < M; i++) {
				if (!GeometryMath::equals(theValues0[i], theValues1[i], theTolerance)) {
					return false;
				}
			}
			return true;
		};

		/// <summary/>
		/// Inspects element-wise equality (using a tolerance value) between two nested std::arrays
		/// 
		/// std::array stores the values instead of using a pointer
		/// Hence, move operation on a std::array has no gain.
		/// The calculations will be performed in a nested loop, in order to optimize the code
		/// </summary>
		template <typename T, size_t M, size_t N>
		static bool equals(
			const std::array<std::array<T, M>, N>& theValues0,
			const std::array<std::array<T, M>, N>& theValues1,
			const double& theTolerance)
		{
			for (int i1 = 0; i1 < N; i1++) {
				for (int i2 = 0; i2 < M; i2++) {
					if (!GeometryMath::equals(theValues0[i1][i2], theValues1[i1][i2], theTolerance)) {
						return false;
					}
				}
			}
			return true;
		};

		/// <summary/>
		/// Inspects equality to zero (using a tolerance value) for all elements in a std::vector
		/// </summary>
		static bool equalsZero(
			const std::vector<double>& theValues,
			const double& theTolerance)
		{
			return std::all_of(
				theValues.cbegin(),
				theValues.cend(),
				[theTolerance](double i) { return std::fabs(i) <= theTolerance; });
		};

		/// <summary/>
		/// Inspects equality to zero (using a tolerance value) for all elements in a std::array
		/// </summary>
		template <typename T, size_t M>
		static bool equalsZero(
			const std::array<T, M>& theValues,
			const double& theTolerance)
		{
			return std::all_of(
				theValues.cbegin(),
				theValues.cend(),
				[theTolerance](double i) { return std::fabs(i) <= theTolerance; });
		};

		/// <summary/>
		/// Inspects equality to zero (using a tolerance value) for all elements in a nested std::array
		/// </summary>
		template <typename T, size_t M, size_t N>
		static bool equalsZero(
			const std::array<std::array<T, M>, N>& theValues,
			const double& theTolerance)
		{
			for (int i = 0; i < N; i++) {
				if (!GeometryMath::equals(theValues[i], 0., theTolerance)) {
					return false;
				}
			}
			return true;
		};

		/// <summary/>
		/// Multiplies each element of the std::array by a factor
		/// </summary>
		template <typename T, size_t M>
		static auto factorizeArray(
			const std::array<T, M>& theValues,
			const double& theFactor)
			-> std::array<T, M>
		{
			std::array<T, M> outValues = { {} };
			std::transform(
				theValues.begin(),
				theValues.end(),
				outValues.begin(),
				[&theFactor](auto& val) {return val * theFactor; });
			return outValues;
		};

		/// <summary/>
		/// Multiplies each element of the nested std::array by a factor
		/// 
		/// std::array stores the values instead of using a pointer
		/// Hence, move operation on a std::array has no gain.
		/// The calculations will be performed in a nested loop, in order to optimize the code
		/// </summary>
		template <typename T, size_t M, size_t N>
		static auto factorizeArray(
			const std::array<std::array<T, M>, N>& theValues,
			const double& theFactor)
			-> std::array<std::array<T, M>, N>
		{
			std::array<std::array<T, M>, N> outValues;
			for (int i1 = 0; i1 < N; i1++) {
				for (int i2 = 0; i2 < M; i2++) {
					outValues[i1][i2] = theValues[i1][i2] * theFactor;
				}
			}
			return outValues;
		};

		/// <summary/>
		/// Performs element-wise summation between two std::arrays
		/// </summary>
		template <typename T, size_t M>
		static auto sumArrays(
			const std::array<T, M>& theValues0,
			const std::array<T, M>& theValues1)
			-> std::array<T, M>
		{
			std::array<T, M> outValues;
			std::transform(
				theValues0.begin(),
				theValues0.end(),
				theValues1.begin(),
				outValues.begin(),
				std::plus<T>());
			return outValues;
		};

		/// <summary/>
		/// Performs element-wise summation between two nested std::arrays
		/// 
		/// std::array stores the values instead of using a pointer
		/// Hence, move operation on a std::array has no gain.
		/// The calculations will be performed in a nested loop, in order to optimize the code
		/// </summary>
		template <typename T, size_t M, size_t N>
		static auto sumArrays(
			const std::array<std::array<T, M>, N>& theValues0,
			const std::array<std::array<T, M>, N>& theValues1)
			-> std::array<std::array<T, M>, N>
		{
			std::array<std::array<T, M>, N> outValues;
			for (int i1 = 0; i1 < N; i1++) {
				for (int i2 = 0; i2 < N; i2++) {
					outValues[i1][i2] = theValues0[i1][i2] + theValues1[i1][i2];
				}
			}
			return outValues;
		};

		/// <summary/>
		/// Performs element-wise subtruction between two std::arrays
		/// </summary>
		template <typename T, size_t M>
		static auto subtructArrays(
			const std::array<T, M>& theValues0,
			const std::array<T, M>& theValues1)
			-> std::array<T, M>
		{
			std::array<T, M> outValues;
			std::transform(
				theValues0.begin(),
				theValues0.end(),
				theValues1.begin(),
				outValues.begin(),
				std::minus<T>());
			return outValues;
		};

		/// <summary/>
		/// Performs element-wise subtruction between two nested std::arrays
		/// 
		/// std::array stores the values instead of using a pointer
		/// Hence, move operation on a std::array has no gain.
		/// The calculations will be performed in a nested loop, in order to optimize the code
		/// </summary>
		template <typename T, size_t M, size_t N>
		static auto subtructArrays(
			const std::array<std::array<T, M>, N>& theValues0,
			const std::array<std::array<T, M>, N>& theValues1)
			-> std::array<std::array<T, M>, N>
		{
			std::array<std::array<T, M>, N> outValues;
			for (int i1 = 0; i1 < N; i1++) {
				for (int i2 = 0; i2 < N; i2++) {
					outValues[i1][i2] = theValues0[i1][i2] - theValues1[i1][i2];
				}
			}
			return outValues;
		};

		/// <summary/>
		/// Converts a 2D nested std::vector with sizes M and N into a nested std:array
		/// </summary>
		template <typename T, size_t M, size_t N>
		static auto convert2DVectorIntoArray(
			const std::vector<std::vector<T, std::allocator<T>>>& theValues)
			-> std::array<std::array<T, M>,N>
		{
			std::array<std::array<T, M>, N> outValues;
			for (int i = 0; i < N; ++i) {
				std::copy(theValues[i].cbegin(), theValues[i].cend(), outValues[i].begin());
			}
			return outValues;
		}

		/// <summary/>
		/// Calculates the multiplication between two MxM matricies both defined by a std::array
		/// </summary>
		template <typename T, size_t M>
		static auto multiplyMatrixToMatrix(
			const std::array<std::array<T, M>, M>& theMatrix0,
			const std::array<std::array<T, M>, M>& theMatrix1)
			-> std::array<std::array<T, M>, M>
		{
			std::array<std::array<T, M>, M> outValues;
			for (int i0 = 0; i0 < M; i0++) {
				for (int i1 = 0; i1 < M; i1++) {
					outValues[i0][i1] = 0.;
					for (int i2 = 0; i2 < M; i2++) {
						outValues[i0][i1] += theMatrix0[i0][i2] * theMatrix1[i2][i1];
					}
				}
			}
			return outValues;
		}

		/// <summary/>
		/// Calculates the multiplication between an MxM matrix and M-sized vector (i.e. 1D array) both defined by a std::array
		/// </summary>
		template <typename T, size_t M>
		static auto multiplyMatrixToVector(
			const std::array<std::array<T, M>, M>& theMatrix,
			const std::array<T, M>& theVector)
			-> std::array<T, M>
		{
			std::array<T, M> outValues = {{}};
			for (int i0 = 0; i0 < M; i0++) {
				outValues[i0] = 0.;
				for (int i1 = 0; i1 < M; i1++) {
					outValues[i0] += theMatrix[i0][i1] * theVector[i1];
				}
			}
			return outValues;
		}

		/// <summary/>
		/// Calculates the determinant of a 3x3 matrix defined by a nested std::array
		/// </summary>
		static double calculateMatrixDeterminantS33(const std::array<std::array<double, 3>, 3>& theMatrix)
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
		static auto calculateMatrixInverseS33(
			const std::array<std::array<double, 3>, 3>& theMatrix,
			const double& theTolerance)
			-> std::array<std::array<double, 3>, 3>
		{
			// Determine determinant
			double determinant = GeometryMath::calculateMatrixDeterminantS33(theMatrix);

			// Check if inversible matrix
			if (std::fabs(determinant) <= theTolerance) throw ZeroDeterminantException();

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

		/// <summary>
		/// Ensures the angle is between -2PI and 2PI
		/// </summary>
		static double normalizeAngle(const double& theAngle)
		{
			double outAngle = theAngle;
			while (std::fabs(outAngle) - M_PI * 2. > TOLERANCE_SENSITIVE) {
				if (outAngle > 0.)
					outAngle -= M_PI * 2.;
				else
					outAngle += M_PI * 2.;
			}
			return outAngle;
		}
	};
}

#endif
