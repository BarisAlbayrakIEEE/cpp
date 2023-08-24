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

	GlobalCoordSystem* GlobalCoordSystem::c_globalCoordSystem;

	/// <summary>
	/// Private constructor of the singleton pattern
	/// </summary>
	GlobalCoordSystem::GlobalCoordSystem()
		: CoordSystem(arrayS3{ 0., 0., 0. }, arrayS3{ 1., 0., 0. }, arrayS3{ 0., 1., 0. }, arrayS3{ 0., 0., 1. })
	{
		c_name = "Global CS";
		c_ID = "Global CS";
		c_isGlobal = true;
	}

	/// <summary>
	/// Static method of the singleton pattern
	/// </summary>
	GlobalCoordSystem* GlobalCoordSystem::getGlobalCoordSystem() {
		if (!c_globalCoordSystem) c_globalCoordSystem = new GlobalCoordSystem();
		return c_globalCoordSystem;
	}

	/// <summary>
	/// The dtor of the singleeton pattern is private to prevent destruction.
	/// </summary>
	GlobalCoordSystem::~GlobalCoordSystem() {
		c_globalCoordSystem = 0;
	}
}
