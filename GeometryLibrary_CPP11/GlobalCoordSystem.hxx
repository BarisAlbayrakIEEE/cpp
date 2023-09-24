/// <summary>
/// GlobalCoordSystem implements the global CS.
/// 
/// All CSs are defined wrt the global CS.
/// 
/// A singleton class
/// 
/// Inherits from CoordSystem.
/// The reason to have a concrete class for the global CS
/// is to have an object simulating the global CS while the program is working.
/// Additionally, ReferenceObjects without a defined CS make use of the global CS.
/// 
/// CAUTION:
///    Despite of being a singleton, the dtor is not removed as:
///       ReferenceObjects has an ownership (via a shared pointer) on CoordSystem objects.
///       ReferenceObject ctors should assign the global CS as the reference CS member
///       if a CoordSystem object is not supplied as an input.
///       However, the GlobalCoordSystem is a singleton
///       that there may exist only one instance and it cannot be deleted
///       which does not allow ownership.
///       There exist two solutions to the problem:
///          1. Remove the destructor as expected for a singleton.
///             Add a static function into CoordSystem class which returns a clone of the GlobalCoordSystem.
///             Call this function from ReferenceObject ctor if the reference CS is not specified.
///             This solution is not efficient as the reference CS is rarely defined in geometrical applications
///             (i.e. mostly work with the global CS)
///             and with this solution,
///             too many CSs are created (cloned) to simulate the global CS.
///          2. Keep the destructor of the GlobalCoordSystem as public default
///             and convert the static GlobalCoordSystem pointer member into a shared pointer.
///             The static member will increment the share count when it is created.
///             Hence, even all reference objects owning the global CS go out of scope,
///             the GlobalCoordSystem object will not be deleted
///             as the static member will still have the ownership on the object.
///             This solution is even better than the static function variable approach
///             (which is the prefered one in C++ singleton design)
///             because, the singleton in our case is a child class and the base class objects are owned by other types.
///       
///       THE 2ND SOLUTION IS SELECTED AS ITS FAR BETTER IN TERMS OF PERFORMANCE
/// 
/// See docstring of CoordSystem.hxx for the CS definition.
/// 
/// See GeometryObject.hxx for the details about this library
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE/cpp.git
/// </summary>

#ifndef _GlobalCoordSystem_HeaderFile
#define _GlobalCoordSystem_HeaderFile

#pragma warning(disable : 4250)

#ifndef _CoordSystem_HeaderFile
#include "CoordSystem.hxx"
#endif
#ifndef _GeometryParameters_HeaderFile
#include "GeometryParameters.hxx"
#endif
#ifndef _GeometryMath_HeaderFile
#include "GeometryMath.hxx"
#endif
#ifndef _GeometryException_HeaderFile
#include "GeometryException.hxx"
#endif
#ifndef _Macros_HeaderFile
#include "Macros.h"
#endif

namespace GeometryNamespace {
	class GlobalCoordSystem : public CoordSystem
	{
		friend class CoordSystem;
		friend class PointBase;
		friend class Point2D;
		friend class Point3D;
		friend class VectorBase;
		friend class Vector20;
		friend class Vector3D;
		friend class Axis;
		friend class Line;
		friend class Circle;
		friend class Plane;

		using CoordSystem::is2D;
		using CoordSystem::is3D;

		// Members
		static std::shared_ptr<GlobalCoordSystem> c_globalCoordSystem;

		// Private ctor of singleton
		GlobalCoordSystem();

	public:
		// copy / move / operators
		GlobalCoordSystem(const GlobalCoordSystem&) = delete;
		GlobalCoordSystem& operator=(const GlobalCoordSystem&) = delete;
		GlobalCoordSystem(GlobalCoordSystem&&) = delete;
		GlobalCoordSystem& operator=(GlobalCoordSystem&&) = delete;
		void setOriginCoords(const std::array<double, 3>& theOriginCoords) = delete;

		// The dtor
		// See the header docstring for details
		~GlobalCoordSystem() final = default;
		static std::shared_ptr<GlobalCoordSystem> getGlobalCoordSystem();
	};
}

#endif
