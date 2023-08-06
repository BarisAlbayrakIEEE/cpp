/// <summary>
/// The abstract base class for the interface inheritance for all types in the library.
/// 
/// CAUTION:
///     Opencascade Technology (OCCT) library shall be embedded to use this library:
///         Download: https://www.opencascade.com/
///         How to: https://www.youtube.com/watch?v=i5zCHArA06E
///     See GeometrySample project in my repository for sample usage of the libraries.
/// 
/// Some pure virtual member functions are not implemented yet (e.g. GetTypeName()).
/// Hence, the definitions for these functions are currently commented out.
/// Later will be implemented.
/// 
/// See GeometryObject.hxx for project definition and main descriptions.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _GeometryAbstractObject_HeaderFile
#define _GeometryAbstractObject_HeaderFile

#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

namespace GeometryNamespace {
	class GeometryAbstractObject
	{
	public:
		Standard_EXPORT virtual void Destroy() = 0;
		//Standard_EXPORT virtual std::string GetTypeName() const = 0;
		//Standard_EXPORT virtual Standard_Size GetSizeOfObject() const = 0;
		//Standard_EXPORT virtual Standard_Size GetDynamicSizeOfObject() const = 0;
	};
}

#endif
