/// <summary>
/// Exceptions defined for the GeometryNamespace
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#ifndef _GeometryException_HeaderFile
#define _GeometryException_HeaderFile

#include <exception>

namespace GeometryNamespace {
	class GeometryException : public std::exception
	{
	public:
		const char* what() const noexcept override
		{
			return "TBD";
		}
	};

	class UncaughtException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Uncaught Exception";
		}
	};

	class AlreadyExistingObjectException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "An object with the input name aalready exists.";
		}
	};

	class ObjectTypeException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "The type of the object does not math the type of the requested object.";
		}
	};

	class NullptrException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Null pointer for an input.";
		}
	};

	class BadInputException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Bad data for an input.";
		}
	};

	class DivideByZeraException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Division by zero";
		}
	};

	class ZeroToleranceException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "The tolerance value set less than or equal to zero which must be positive real value.";
		}
	};

	class ZeroDeterminantException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Determinant of the matrix is zerom. Inverse matrix does not exist.";
		}
	};

	class ArraySizeException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Array size is zero ar not same as the expected";
		}
	};

	class MaxIterationReachedException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "No solution. Max iteration has reached.";
		}
	};

	class DimensionalityException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "The dimensionality of the instance is not suitable for the requested operation.";
		}
	};

	class ZeroDimensionException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "The instance is defined in 2D space.";
		}
	};

	class GeometricalMismatchException : public GeometryException
	{
	public:
		const char* what() const noexcept override
		{
			return "Two or more entities in the operation do not match (e.g. a point does not lie on the plane or non-perpandicular vectors for a CS).";
		}
	};

	class CoordSystemMismatchException : public GeometricalMismatchException
	{
	public:
		const char* what() const noexcept final
		{
			return "Reference CSs are not the same.";
		}
	};

	class NoIntersectionException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "No intersection between the two components.";
		}
	};

	class ZeroVectorException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Vector magnitude is zero. nullptr vector";
		}
	};

	class CoincidenceException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Coincident items";
		}
	};

	class ColinearPointsException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Points are colinear.";
		}
	};

	class ParallelAxisException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "Axes are parallel.";
		}
	};

	class AssymptoticLineException : public GeometryException
	{
	public:
		const char* what() const noexcept final
		{
			return "The line is parallel to one of the axes.";
		}
	};
}

#endif
