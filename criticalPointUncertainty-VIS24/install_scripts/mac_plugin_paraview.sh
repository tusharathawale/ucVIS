#!/bin/bash
set -e

build_jobs=6
mkdir -p mac
cd mac

HERE=`pwd`
source $HERE/../settings.sh

SOFTWARE_SRC_DIR="$HERE/src"
SOFTWARE_BUILD_DIR="$HERE/build"
SOFTWARE_INSTALL_DIR="$HERE/install"

mkdir -p $SOFTWARE_SRC_DIR
mkdir -p $SOFTWARE_BUILD_DIR
mkdir -p $SOFTWARE_INSTALL_DIR


echo "====> Installing vtk-m"
VTKM_SRC_DIR="$SOFTWARE_SRC_DIR/vtk-m"
VTKM_BUILD_DIR="$SOFTWARE_BUILD_DIR/vtk-m"
VTKM_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/vtk-m"

# check the install dir
if [ -d $VTKM_INSTALL_DIR ]; then
    echo "====> skip, $VTKM_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else
    echo $VTKM_SRC_DIR
    echo $VTKM_BUILD_DIR
    echo $VTKM_INSTALL_DIR
    # check vktm source dir
    if [ ! -d $VTKM_SRC_DIR ]; then
    # clone the source
    cd $SOFTWARE_SRC_DIR
    git clone $VTKM_REPO
    cd $VTKM_SRC_DIR
    #git checkout v2.0.0-rc1
    #git checkout $VTKM_VERSION    
fi
    
    cd $HERE

    # build and install
    echo "**** Building vtk-m"

    # TODO, the gpu version can be different here
    # we only use the cpu version here

    cmake -B ${VTKM_BUILD_DIR} -S ${VTKM_SRC_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DVTKm_USE_DEFAULT_TYPES_FOR_ASCENT=ON \
    -DVTKm_USE_DOUBLE_PRECISION=ON \
    -DVTKm_USE_64BIT_IDS=OFF \
    -DCMAKE_INSTALL_PREFIX=${VTKM_INSTALL_DIR} \
    -DVTKm_ENABLE_MPI=OFF \
    -DVTKm_ENABLE_TBB=OFF \
    -DVTKm_ENABLE_LOGGING=ON \
    -DVTKm_ENABLE_TESTING=OFF \
    -DVTKm_ENABLE_RENDERING=OFF 

    cmake --build ${VTKM_BUILD_DIR} -j${build_jobs}

    echo "**** Installing vtk-m"
    cmake --install ${VTKM_BUILD_DIR}
fi

echo "====> Installing vtk-m, ok"


echo "====> Installing EasyLinalg"
EASY_LINALG_SRC_DIR="$SOFTWARE_SRC_DIR/EasyLinalg"
EASY_LINALG_INSTALL_DIR="$HERE/../../ucvworklet/linalg/EasyLinalg/"

rm -rf $EASY_LINALG_SRC_DIR
cd $SOFTWARE_SRC_DIR
git clone $EASY_LINALG_REPO

# move include dir to correct place

# clean old dir if it exist
if [ -d $EASY_LINALG_INSTALL_DIR ]; then
    rm -rf $EASY_LINALG_INSTALL_DIR
fi

mkdir -p $EASY_LINALG_INSTALL_DIR

# move files to new dir
cp EasyLinalg/StaticMemTemplate/include/* $EASY_LINALG_INSTALL_DIR
# clean source files
rm -rf $EASY_LINALG_SRC_DIR

echo "====> Installing EasyLinalg, ok"


echo "====> build UCV"
# the only have build dir without the install dir
# the install dir is same with the build dir
UCV_SRC_DIR=$HERE/../../
# use the install dir as the build dir
UCV_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/UCV"
rm -rf $UCV_INSTALL_DIR

#build the latest paraview (paraview 1.2 for this)
#if [ -d $UCV_INSTALL_DIR ]; then
#    echo "====> skip, $UCV_INSTALL_DIR already exists," \
#             "please remove it if you want to reinstall it"
#else
    
    #-DCMAKE_BUILD_TYPE=Release \
    # using Debug mode for enabling the assert option in testing
    # comment out paraview dir if not install plugin
    # if using paraview plugin build, just comment out the vtkm dir
    cmake -B ${UCV_INSTALL_DIR} -S ${UCV_SRC_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DVTKm_DIR=/Users/tm8/tusharathawale/projects/software/criticalPointUncertainty/UCV/install_scripts/mac/install/paraview/lib/cmake/paraview-5.13/vtk/vtkm \
    -DParaView_DIR=/Users/tm8/tusharathawale/projects/software/criticalPointUncertainty/UCV/install_scripts/mac/install/paraview/lib/cmake/paraview-5.13 \
    -DQt5_DIR=/usr/local/Cellar/qt@5/5.15.13_1/lib/cmake/Qt5

    cd $HERE

    # build and install
    echo "**** Building UCV"
    cmake --build ${UCV_INSTALL_DIR} -j${build_jobs}
#fi

# not sure why the libvtkmdiympi.so is not included during the build process
echo "try to add library path by executing:"
echo "export LD_LIBRARY_PATH=${VTKM_INSTALL_DIR}/lib:\${LD_LIBRARY_PATH}"
