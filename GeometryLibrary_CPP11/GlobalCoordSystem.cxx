// baris.albayrak.ieee@gmail.com

#include "GeometryObject.hxx"
#include "ReferenceObject.hxx"
#include "CoordSystem.hxx"
#include "GlobalCoordSystem.hxx"
#include "PointBase.hxx"
#include "Point2D.hxx"
#include "Point3D.hxx"
#include "VectorBase.hxx"
#include "Vector2D.hxx"
#include "Vector3D.hxx"
#include "Axis.hxx"
#include "Line.hxx"
#include "Circle.hxx"
#include "Plane.hxx"

namespace GeometryNamespace {

	std::shared_ptr<GlobalCoordSystem> GlobalCoordSystem::c_globalCoordSystem;

	/// <summary>
	/// Private constructor of the singleton pattern
	/// </summary>
	GlobalCoordSystem::GlobalCoordSystem()
		: CoordSystem(
			std::array<double, 3>{{ 0., 0., 0. }},
			std::array<double, 3>{{ 1., 0., 0. }},
			std::array<double, 3>{{ 0., 1., 0. }},
			std::array<double, 3>{{ 0., 0., 1. }})
	{
		setName("Global CS");
		setID(0);
		setIsGlobal(true);
	}

	/// <summary>
	/// Static method of the singleton pattern
	/// </summary>
	std::shared_ptr<GlobalCoordSystem> GlobalCoordSystem::getGlobalCoordSystem() {
		if (!c_globalCoordSystem) c_globalCoordSystem = std::shared_ptr<GlobalCoordSystem>(new GlobalCoordSystem());
		return c_globalCoordSystem;
	}
}
