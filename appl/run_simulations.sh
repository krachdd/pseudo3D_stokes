#!/bin/sh
# Define your module path here
MODULEPATH="$(dirname "$(pwd)")"
echo "Modulepath: $MODULEPATH"
DUMUXPATH="${MODULEPATH}/build-cmake/appl"
EXECUTABLE="pseudo3D_stokes"
export PYTHONPATH=$PYTHONPATH:${MODULEPATH}/preprocessing/localdrag/

# TODO list the relevant input params here, with explanation

# runs 
python3 runStokesGeneric.py -dir test_singlePrecipitate_2d_total -vs 1e-6 -height 36.0e-6 -dumux $DUMUXPATH -exec $EXECUTABLE -get_phi -vtu -govEq brinkman
