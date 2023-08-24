/// <summary>
/// GlobalCoordSystem is a singleton class that implements the global CS.
/// 
/// All CSs are defined wrt the global CS.
/// 
/// Inherits from CoordSystem.
/// The reason to have a concrete class for the global CS
/// is to emphasize the importance of the global CS for the users of the library.
/// Additionally, the default ctors make use of the global CS.
/// 
/// See docstring of CoordSystem.hxx for the CS definition.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
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
	class GeometryObject;
	class ReferenceObject;
	class CoordSystem;
	class GlobalCoordSystem;
	class PointBase;
	class Point2D;
	class Point3D;
	class VectorBase;
	class Vector2D;
	class Vector3D;
	class Axis;
	class Line;
	class Circle;
	class Plane;

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

	private:
		GlobalCoordSystem();

	private:
		// Members
		static GlobalCoordSystem* c_globalCoordSystem;

	public:
		// copy / move / operators
		GlobalCoordSystem(const GlobalCoordSystem&) = delete;
		void operator=(const GlobalCoordSystem&) = delete;
		GlobalCoordSystem(GlobalCoordSystem&&) = delete;
		void operator=(GlobalCoordSystem&&) = delete;

	private:
		~GlobalCoordSystem();

	public:
		static GlobalCoordSystem* getGlobalCoordSystem();

	};
}

#endif
