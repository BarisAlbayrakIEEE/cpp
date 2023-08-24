/// <summary>
/// Includes some simple math functions.
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#ifndef _GeometryMath_hxx_
#define _GeometryMath_hxx_

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

namespace GeometryNamespace {
	struct GeometryMath
	{
	public:
		static double interpolate(const double& theValueO, const double& theValue1, const double& theFactor);
		static bool equals(const double& theValueO, const double& theValue1, const double& theTolerance);
		template <typename T, size_t M>
		static bool equals(
			std::array<T, M> theValues0,
			std::array<T, M> theValues1,
			const double& theTolerance);
		template <typename T, size_t M, size_t N>
		static bool equals(
			std::array<std::array<T, M>, N> theValues0,
			std::array<std::array<T, M>, N> theValues1,
			const double& theTolerance);
		static bool equalsZero(
			std::vector<double> theValues,
			const double& theTolerance);
		template <typename T, size_t M>
		static bool equalsZero(
			std::array<T, M> theValues,
			const double& theTolerance);
		template <typename T, size_t M, size_t N>
		static bool equalsZero(
			std::array<std::array<T, M>, N> theValues,
			const double& theTolerance);
		template <typename T, size_t M, size_t N>
		static std::array<std::array<T, M>, N> copyNestedArray(std::array<std::array<T, M>, N> theValues);
		template <typename T, size_t M>
		static std::array<T, M> factorizeArray(
			std::array<T, M> theValues,
			const double& theFactor);
		template <typename T, size_t M, size_t N>
		static std::array<std::array<T, M>, N> factorizeArray(
			std::array<std::array<T, M>, N> theValues,
			const double& theFactor);
		template <typename T, size_t M>
		static std::array<T, M> sumArrays(
			std::array<T, M> theValues0,
			std::array<T, M> theValues1);
		template <typename T, size_t M, size_t N>
		static std::array<std::array<T, M>, N> sumArrays(
			std::array<std::array<T, M>, N> theValues0,
			std::array<std::array<T, M>, N> theValues1);
		template <typename T, size_t M>
		static std::array<T, M> subtructArrays(
			std::array<T, M> theValues0,
			std::array<T, M> theValues1);
		template <typename T, size_t M, size_t N>
		static std::array<std::array<T, M>, N> subtructArrays(
			std::array<std::array<T, M>, N> theValues0,
			std::array<std::array<T, M>, N> theValues1);

		static arrayS3 convertArrayS2ToS3(const arrayS2& theValues);
		static arrayS2 convertVectorToArray1DS2(const vectorInput1D& theValues);
		static arrayS3 convertVectorToArray1DS3(const vectorInput1D& theValues);
		static arrayS4 convertVectorToArray1DS4(const vectorInput1D& theValues);
		static arrayS32 convertVectorToArray2DS32(const vectorInput2D& theValues);
		static arrayS33 convertVectorToArray2DS33(const vectorInput2D& theValues);
		static arrayS33 multiplyMatrixToMatrixS3(const arrayS33& theValues0, const arrayS33& theValues1);
		static arrayS3 multiplyMatrixToVectorS33S3(const arrayS33& theValues0, const arrayS3& theValues1);
		static double calculateMatrixDeterminantS33(const arrayS33& theMatrix);
		static arrayS33 calculateMatrixInverseS33(const arrayS33& theMatrix, const double& theTolerance);
		static double normalizeAngle(const double& theAngle);
	};
};

#endif
