/// <summary>
/// GlobalCoordSystem is a singleton class that implements the global CS.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// All CSs are defined wrt the global CS.
/// 
/// Inherits from CoordSystem.
/// The reason to have a concrete class for the global CS
/// is to emphasize the importaance of the global CS for the users of the library.
/// 
/// See docstring of CoordSystem.hxx for the CS definition.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#ifndef _GlobalCoordSystem_HeaderFile
#define _GlobalCoordSystem_HeaderFile

#pragma warning(disable : 4250)

#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Handle_HeaderFile
#include <Standard_Handle.hxx>
#endif
#ifndef _Standard_Type_HeaderFile
#include <Standard_Type.hxx>
#endif
#ifndef _Standard_Size_HeaderFile
#include <Standard_Size.hxx>
#endif
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

	DEFINE_STANDARD_HANDLE(GlobalCoordSystem, CoordSystem)

	class GlobalCoordSystem : public virtual GeometryAbstractObject, public CoordSystem
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

	public:
		/*
		 * DEF: Modifications on OCCT
		*/
		//void* operator new(size_t theSize, void* theAdress) { return theAdress; }
		//void* operator new(size_t theSize) { return Standard::Allocate(theSize); }
		//void operator delete(void* theAdress) { if (theAdress) Standard::Free((Standard_Address&)theAdress); }
		Standard_EXPORT virtual Standard_Size getSizeOfObject() const {
			return sizeof(GlobalCoordSystem);
		}
		/*
		 * DEF END: Modifications on OCCT
		*/

	public:
		Standard_EXPORT void Destroy();

	private:
		Standard_EXPORT GlobalCoordSystem();

	private:
		// Members
		static GlobalCoordSystem* c_globalCoordSystem;

	public:
		// copy / move / operators
		Standard_EXPORT GlobalCoordSystem(const GlobalCoordSystem&) = delete;
		Standard_EXPORT void operator=(const GlobalCoordSystem&) = delete;
		Standard_EXPORT GlobalCoordSystem(GlobalCoordSystem&&) = delete;
		Standard_EXPORT void operator=(GlobalCoordSystem&&) = delete;

	private:
		Standard_EXPORT ~GlobalCoordSystem();

	public:
		Standard_EXPORT static GlobalCoordSystem* getGlobalCoordSystem();

	public:
		DEFINE_STANDARD_RTTIEXT(GlobalCoordSystem, CoordSystem)
	};
}

#endif
