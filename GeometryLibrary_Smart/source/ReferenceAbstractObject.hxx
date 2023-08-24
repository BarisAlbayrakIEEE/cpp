/// <summary>
/// The abstract base class for the interface inheritance
/// for the reference type objects (i.e. Point and Vector).
/// 
/// Some pure virtual member functions are not implemented yet.
/// Later will be implemented.
/// 
/// author: baris.albayrak.ieee@gmail.com
/// github: https://github.com/BarisAlbayrakIEEE?tab=repositories
/// </summary>

#pragma warning(disable : 4290)

#ifndef _ReferenceAbstractObject_HeaderFile
#define _ReferenceAbstractObject_HeaderFile

namespace GeometryNamespace {
	class ReferenceAbstractObject
	{
	public:
		virtual bool is2D() const = 0;
		virtual bool is3D() const = 0;
	};
}

#endif
