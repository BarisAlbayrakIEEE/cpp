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

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>

#ifndef _GeometryException_HeaderFile
#include "GeometryException.hxx"
#endif
#ifndef _GeometryParameters_HeaderFile
#include "GeometryParameters.hxx"
#endif

namespace GeometryNamespace {
	class GeometryMath
	{
	public:
		// Comparison functions
		template<typename T>
		static bool zero_g(T&& theValue)
		{
			return std::fabs(std::forward<T>(theValue)) < GeometryParameters::getToleranceGeneral();
		};
		template<typename T>
		static bool zero_s(T&& theValue)
		{
			return std::fabs(std::forward<T>(theValue)) < GeometryParameters::getToleranceSensitive();
		};
		template<typename T>
		static bool GTEZ_g(T&& theValue)
		{
			return std::forward<T>(theValue) >= -GeometryParameters::getToleranceGeneral();
		};
		template<typename T>
		static bool GTEZ_s(T&& theValue)
		{
			return std::forward<T>(theValue) >= -GeometryParameters::getToleranceSensitive();
		};
		template<typename T>
		static bool LTEZ_g(T&& theValue)
		{
			return std::forward<T>(theValue) <= GeometryParameters::getToleranceGeneral();
		};
		template<typename T>
		static bool LTEZ_s(T&& theValue)
		{
			return std::forward<T>(theValue) <= GeometryParameters::getToleranceSensitive();
		};
		template<typename T>
		static bool equal_g(T&& theValue1, T&& theValue2)
		{
			return zero_g(std::forward<T>(theValue1) - std::forward<T>(theValue2));
		};
		template<typename T>
		static bool equal_s(T&& theValue1, T&& theValue2)
		{
			return zero_s(std::forward<T>(theValue1) - std::forward<T>(theValue2));
		};
		template<typename T>
		static bool absoluteEqual_g(T&& theValue1, T&& theValue2)
		{
			return zero_g(std::fabs(std::forward<T>(theValue1)) - std::fabs(std::forward<T>(theValue2)));
		};
		template<typename T>
		static bool absoluteEqual_s(T&& theValue1, T&& theValue2)
		{
			return zero_s(std::fabs(std::forward<T>(theValue1)) - std::fabs(std::forward<T>(theValue2)));
		};
		template <typename T, typename Iterator>
		static bool allAbsoluteEqual_g(Iterator theBegin0, Iterator theEnd0, Iterator theBegin1)
		{
			return std::equal(
				theBegin0,
				theEnd0,
				theBegin1,
				[](T i, T j) { return absoluteEqual_g(i, j); });
		};
		template <typename T, typename Iterator>
		static bool allAbsoluteEqual_s(Iterator theBegin0, Iterator theEnd0, Iterator theBegin1)
		{
			return std::equal(
				theBegin0,
				theEnd0,
				theBegin1,
				[](T i, T j) { return absoluteEqual_s(i, j); });
		};

		// Basic math functions
		template<typename T>
		static auto interpolate(const T& theValue1, const T& theValue2, const T& theFactor) -> T
		{
			return theValue1 + (theValue2 - theValue1) * theFactor;
		};

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
			std::array<T, M> outValues = { {} };
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
		static auto calculateMatrixDeterminantS33(const std::array<std::array<double, 3>, 3>& theMatrix) -> double;

		/// <summary/>
		/// Calculates the inverse of a 3x3 matrix defined by a nested std::array
		/// </summary>
		static auto calculateMatrixInverseS33(
			const std::array<std::array<double, 3>, 3>& theMatrix)
			-> std::array<std::array<double, 3>, 3>;
	};
}

#endif
