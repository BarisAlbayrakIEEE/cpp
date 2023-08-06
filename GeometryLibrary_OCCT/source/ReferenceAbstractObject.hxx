/// <summary>
/// The abstract base class for the interface inheritance
/// for the reference type objects (i.e. Point and Vector).
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// Some pure virtual member functions are not implemented yet.
/// Later will be implemented.
/// 
/// See GeometryAbstractObject.hxx for the abstract base class of the library.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _ReferenceAbstractObject_HeaderFile
#define _ReferenceAbstractObject_HeaderFile

#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

namespace GeometryNamespace {
	class ReferenceAbstractObject : public GeometryAbstractObject
	{
	public:
		Standard_EXPORT virtual bool is2D() const = 0;
		Standard_EXPORT virtual bool is3D() const = 0;
	};
}

#endif
