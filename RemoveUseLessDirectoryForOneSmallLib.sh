#!/bin/sh

#rm -rf .git
rm -rf Accelerators Charts Deprecated Domains Examples Geovis GUISupport Imaging Infovis Interaction Rendering Testing Views Web Wrapping
rm -rf Common/Color Common/ComputationalGeometry
rm -rf Filters/A* Filters/F* Filters/G* Filters/H* Filters/I*
rm -rf Filters/M* Filters/P* Filters/R* Filters/S* Filters/T* Filters/V*
rm -rf IO/A* IO/E* IO/F* IO/GDAL IO/GeoJSON IO/I* IO/LSD* IO/M* IO/N* IO/O* IO/P* IO/SQL IO/Q* IO/V* IO/Xd*
rm -rf ThirdParty/a* ThirdParty/A* ThirdParty/exodus* ThirdParty/f* ThirdParty/g* ThirdParty/j* ThirdParty/l*
rm -rf ThirdParty/m* ThirdParty/n* ThirdParty/o* ThirdParty/p* ThirdParty/s* ThirdParty/SixPython ThirdParty/T* ThirdParty/t*
rm -rf hirdParty/hdf5 ThirdParty/verdict ThirdParty/VPIC ThirdParty/xdmf2 ThirdParty/xdmf3 ThirdParty/ZopeInterface
rm -rf Utilities/Benchmarks Utilities/D* Utilities/E* Utilities/G* Utilities/M* Utilities/oct* Utilities/P* Utilities/R* Utilities/S* Utilities/U*  Utilities/vtkTcl*
rm -rf .ExternalData

find . -name Testing | xargs rm -rf
