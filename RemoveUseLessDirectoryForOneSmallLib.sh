#!/bin/bash

shopt -s extglob

#rm -rf .git
rm -rf .ExternalData Accelerators Charts Deprecated Documentation Domains Examples Geovis GUISupport Imaging Infovis Interaction Rendering Testing Views Web Wrapping
cd Common
rm -rf !(ComputationalGeometry|Core|DataModel|ExecutionModel|Math|Misc|System|Transforms)
cd ../Filters
rm -rf !(Core|Extraction|General|Geometry|Sources)
cd ../IO
rm -rf !(Core|Geometry|Legacy|ParallelXML|XML|XMLParser)
cd ../Parallel
rm -rf !(Core|MPI|MPI4Py)
cd ../ThirdParty
rm -rf !(expat|utf8|zlib|update-common.sh)
cd ../Utilities
rm -rf !(HashSource|KWSys|KWIML|OutputWindowProcess|*in|*py)
cd ..

find . -name Testing | xargs rm -rf

echo "Next step to do manually:"
echo "cd /path/ShaPo/thirdparty/VTK"
echo "rsync -av --exclude=README.shapo --exclude=VERSION --delete ~/dev/VTK/ ./ | grep deleting | cut -b10- | xargs git rm"
echo "Then: fix potential module dependencies manually & remove some filters"
