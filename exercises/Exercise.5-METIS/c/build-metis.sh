#!/bin/bash

# ======================= #
# project directory paths
# ======================= #
PROJECT_ROOT="$(pwd)"
INSTALL_DIRECTORY=${PROJECT_ROOT}/install
INSTALL_METIS_DIRECTORY=${INSTALL_DIRECTORY}/metis

BUILD_TYPE="Release"

cd metis-5.0.3
rm -rf buildW
rm -rf ${INSTALL_METIS_DIRECTORY}
mkdir buildW

pushd buildW
cmake -D CMAKE_INSTALL_PREFIX:PATH=${PROJECT_ROOT}/metis-5.0.3/install \
      -D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}                         \
      -D GKLIB_PATH:PATH=${PROJECT_ROOT}/metis-5.0.3/GKlib             \
      -D SHARED:BOOL=TRUE                                              \
      -G "Unix Makefiles" ../

make -j4 install
popd

