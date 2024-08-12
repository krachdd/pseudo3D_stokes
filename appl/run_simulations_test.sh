#!/bin/sh
# Define your module path here
MODULEPATH="$(dirname "$(dirname "$(pwd)")")"
echo $MODULEPATH
DUMUXPATH="${MODULEPATH}/build-cmake/appl"
EXECUTABLE="pseudo3D_stokes"
export PYTHONPATH=$PYTHONPATH:${MODULEPATH}/preprocessing/localdrag/

# runs 
python3 runStokesGeneric.py -dir test_singlePrecipitate_2d_total -vs 1e-6 -height 36.0e-6 -dumux $DUMUXPATH -exec $EXECUTABLE -get_phi -vtu
