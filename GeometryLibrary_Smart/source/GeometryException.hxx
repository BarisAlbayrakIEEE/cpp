/// <summary>
/// Exceptions defined for the GeometryNamespace
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#ifndef _GeometryException_hxx_
#define _GeometryException_hxx_

#include <iostream>
using namespace std;

namespace GeometryNamespace {
	class GeometryException : public std::exception
	{
	public:
		const char* what() const noexcept
		{
			return "TBD";
		}
	};

	class UncaughtException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Uncaught Exception";
		}
	};

	class AlreadyExistingObjectException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "An object with the input name aalready exists.";
		}
	};

	class ObjectTypeException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "The type of the object does not math the type of the requested object.";
		}
	};

	class NullptrException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Null pointer for an input.";
		}
	};

	class BadInputException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Bad data for an input.";
		}
	};

	class DivideByZeraException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Division by zero";
		}
	};

	class ZeroToleranceException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "The tolerance value set less than or equal to zero which must be positive real value.";
		}
	};

	class ZeroDeterminantException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Determinant of the matrix is zerom. Inverse matrix does not exist.";
		}
	};

	class ArraySizeException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Array size is zero ar not same as the expected";
		}
	};

	class MaxIterationReachedException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "No solution. Max iteration has reached.";
		}
	};

	class DimensionalityException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "The dimensionality of the instance is not suitable for the requested operation.";
		}
	};

	class ZeroDimensionException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "The instance is defined in 2D space.";
		}
	};

	class GeometricalMismatchException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Two or more entities in the operation do not match (e.g. a point does not lie on the plane or non-perpandicular vectors for a CS).";
		}
	};

	class CoordSystemMismatchException : public GeometricalMismatchException
	{
	public:
		const char* what() const noexcept
		{
			return "Reference CSs are not the same.";
		}
	};

	class NoIntersectionException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "No intersection between the two components.";
		}
	};

	class ZeroVectorException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Vector magnitude is zero. nullptr vector";
		}
	};

	class CoincidenceException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Coincident items";
		}
	};

	class ColinearPointsException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Points are colinear.";
		}
	};

	class ParallelAxisException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "Axes are parallel.";
		}
	};

	class AssymptoticLineException : public GeometryException
	{
	public:
		const char* what() const noexcept
		{
			return "The line is parallel to one of the axes.";
		}
	};
}

#endif
