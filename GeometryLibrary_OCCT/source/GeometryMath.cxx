// baris.albayrak.ieee@gmail.com

#include "GeometryMath.hxx"

namespace GeometryNamespace {
	/// <summary>
	/// theFactor: Betwen 0. and 1.
	///	Defines how much of the distance between the two values is to be travelled
	/// </summary>
	double GeometryMath::interpolate(const double& theValue0, const double& theValue1, const double& theFactor)
	{
		return theValue0 + (theValue1 - theValue0) * theFactor;
	}

	bool GeometryMath::equals(const double& theValue0, const double& theValue1, const double& theTolerance)
	{
		if (std::fabs(theValue0 - theValue1) < theTolerance) return true;
		return false;
	}

	template <typename T, size_t M>
	bool GeometryMath::equals(std::array<T, M> theValues0, std::array<T, M> theValues1, const double& theTolerance) {
		for (int i = 0; i < M; i++) {
			if (!GeometryMath::equals(theValues0[i], theValues1[i], theTolerance)) {
				return false;
			}
		}
		return true;
	};

	template <typename T, size_t M, size_t N>
	bool GeometryMath::equals(
		std::array<std::array<T, M>, N> theValues0,
		std::array<std::array<T, M>, N> theValues1,
		const double& theTolerance) {
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 < M; i2++) {
				if (!GeometryMath::equals(theValues0[i1][i2], theValues1[i1][i2], theTolerance)) {
					return false;
				}
			}
		}
		return true;
	};

	template <typename T, size_t M>
	bool GeometryMath::equalsZero(std::array<T, M> theValues, const double& theTolerance) {
		for (int i = 0; i < M; i++) {
			if (!GeometryMath::equals(theValues[i], 0., theTolerance)) {
				return false;
			}
		}
		return true;
	};

	template <typename T, size_t M, size_t N>
	bool GeometryMath::equalsZero(std::array<std::array<T, M>, N> theValues, const double& theTolerance) {
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 < M; i2++) {
				if (!GeometryMath::equals(theValues[i1][i2], 0., theTolerance)) {
					return false;
				}
			}
		}
		return true;
	};

	template <typename T, size_t M, size_t N>
	std::array<std::array<T, M>, N> GeometryMath::copyNestedArray(std::array<std::array<T, M>, N> theValues) {
		std::array<std::array<T, M>, N> outValues;
		for (int i = 0; i < N; ++i) {
			std::copy(std::begin(theValues[i]), std::end(theValues[i]), std::begin(outValues[i]));
		}
		return outValues;
	};

	template <typename T, size_t M>
	std::array<T, M> GeometryMath::factorizeArray(
		std::array<T, M> theValues,
		const double& theFactor) {
		std::array<T, M> outValues;
		for (int i = 0; i < M; i++) {
			outValues[i] = theFactor * theValues[i];
		}
		return outValues;
	};

	template <typename T, size_t M, size_t N>
	std::array<std::array<T, M>, N> GeometryMath::factorizeArray(
		std::array<std::array<T, M>, N> theValues,
		const double& theFactor) {
		std::array<std::array<T, M>, N> outValues;
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 < M; i2++) {
				outValues[i1][i2] = theFactor * theValues[i1][i2];
			}
		}
		return outValues;
	};

	template <typename T, size_t M>
	std::array<T, M> GeometryMath::sumArrays(
		std::array<T, M> theValues0,
		std::array<T, M> theValues1) {
		std::array<T, M> outValues;
		for (int i = 0; i < M; i++) {
			outValues[i] = theValues0[i] + theValues1[i];
		}
		return outValues;
	};

	template <typename T, size_t M, size_t N>
	std::array<std::array<T, M>, N> GeometryMath::sumArrays(
		std::array<std::array<T, M>, N> theValues0,
		std::array<std::array<T, M>, N> theValues1) {
		std::array<std::array<T, M>, N> outValues;
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 < M; i2++) {
				outValues[i1][i2] = theValues0[i1][i2] + theValues1[i1][i2];
			}
		}
		return true;
	};

	template <typename T, size_t M>
	std::array<T, M> GeometryMath::subtructArrays(
		std::array<T, M> theValues0,
		std::array<T, M> theValues1) {
		std::array<T, M> outValues;
		for (int i = 0; i < M; i++) {
			outValues[i] = theValues0[i] - theValues1[i];
		}
		return outValues;
	};

	template <typename T, size_t M, size_t N>
	std::array<std::array<T, M>, N> GeometryMath::subtructArrays(
		std::array<std::array<T, M>, N> theValues0,
		std::array<std::array<T, M>, N> theValues1) {
		std::array<std::array<T, M>, N> outValues;
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 < M; i2++) {
				outValues[i1][i2] = theValues0[i1][i2] - theValues1[i1][i2];
			}
		}
		return true;
	};

	arrayS3 GeometryMath::convertArrayS2ToS3(const arrayS2& theValues) {
		arrayS3 s3{ arrayS3{} };
		for (int i = 0; i < DIMENSIONS::D2; ++i)
		{
			s3[i] = theValues[i];
		}
		return s3;
	}

	arrayS2 GeometryMath::convertVectorToArray1DS2(const vectorInput1D& theValues) {
		if (theValues.size() != 2) throw ArraySizeException();

		arrayS2 outValues{ arrayS2{} };
		for (int i = 0; i < 2; ++i) {
			outValues[i] = theValues[i];
		}
		return outValues;
	}

	arrayS3 GeometryMath::convertVectorToArray1DS3(const vectorInput1D& theValues) {
		if (theValues.size() != 3) throw ArraySizeException();

		arrayS3 outValues{ arrayS3{} };
		for (int i = 0; i < 3; ++i) {
			outValues[i] = theValues[i];
		}
		return outValues;
	}

	arrayS4 GeometryMath::convertVectorToArray1DS4(const vectorInput1D& theValues) {
		if (theValues.size() != 4) throw ArraySizeException();

		arrayS4 outValues{ arrayS4{} };
		for (int i = 0; i < 4; ++i) {
			outValues[i] = theValues[i];
		}
		return outValues;
	}

	arrayS32 GeometryMath::convertVectorToArray2DS32(const vectorInput2D& theValues) {
		if (theValues.size() != 3) throw ArraySizeException();
		if (theValues[0].size() != 2) throw ArraySizeException();
		if (theValues[1].size() != 2) throw ArraySizeException();
		if (theValues[2].size() != 2) throw ArraySizeException();

		arrayS32 outValues{ arrayS32{} };
		for (int i0 = 0; i0 < 3; ++i0) {
			for (int i1 = 0; i1 < 2; ++i1) {
				outValues[i0][i1] = theValues[i0][i1];
			}
		}
		return outValues;
	}

	arrayS33 GeometryMath::convertVectorToArray2DS33(const vectorInput2D& theValues) {
		if (theValues.size() != 3) throw ArraySizeException();
		if (theValues[0].size() != 3) throw ArraySizeException();
		if (theValues[1].size() != 3) throw ArraySizeException();
		if (theValues[2].size() != 3) throw ArraySizeException();

		arrayS33 outValues{ arrayS33{} };
		for (int i0 = 0; i0 < 3; ++i0) {
			for (int i1 = 0; i1 < 3; ++i1) {
				outValues[i0][i1] = theValues[i0][i1];
			}
		}
		return outValues;
	}

	arrayS33 GeometryMath::multiplyMatrixToMatrixS3(const arrayS33& theMatrix0, const arrayS33& theMatrix1) {
		arrayS33 outValues{ arrayS33{} };
		for (int i0 = 0; i0 < theMatrix0.size(); i0++) {
			for (int i1 = 0; i1 < theMatrix0.size(); i1++) {
				outValues[i0][i1] = 0.;
				for (int i2 = 0; i2 < theMatrix0.size(); i2++) {
					outValues[i0][i1] += theMatrix0[i0][i2] * theMatrix1[i2][i1];
				}
			}
		}
		return outValues;
	}

	arrayS3 GeometryMath::multiplyMatrixToVectorS33S3(const arrayS33& theMatrix, const arrayS3& theVector)
	{
		arrayS3 outValues{ arrayS3{} };
		for (int i0 = 0; i0 < theMatrix.size(); i0++) {
			outValues[i0] = 0.;
			for (int i1 = 0; i1 < theMatrix.size(); i1++) {
				outValues[i0] += theMatrix[i0][i1] * theVector[i1];
			}
		}
		return outValues;
	}

	double GeometryMath::calculateMatrixDeterminantS33(const arrayS33& theMatrix)
	{
		double outDeterminant = 0.;
		outDeterminant += theMatrix[0][0] * (theMatrix[1][1] * theMatrix[2][2] - theMatrix[1][2] * theMatrix[2][1]);
		outDeterminant -= theMatrix[0][1] * (theMatrix[1][0] * theMatrix[2][2] - theMatrix[1][2] * theMatrix[2][0]);
		outDeterminant += theMatrix[0][2] * (theMatrix[1][0] * theMatrix[2][1] - theMatrix[1][1] * theMatrix[2][0]);
		return outDeterminant;
	}

	arrayS33 GeometryMath::calculateMatrixInverseS33(const arrayS33& theMatrix, const double& theTolerance)
	{
		// Determine determinant
		double determinant = GeometryMath::calculateMatrixDeterminantS33(theMatrix);

		// Check if inversible matrix
		if (std::fabs(determinant) <= theTolerance) throw ZeroDeterminantException();

		// The inverse matrix
		arrayS33 outInverseMatrix;
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
	double GeometryMath::normalizeAngle(const double& theAngle)
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
}
