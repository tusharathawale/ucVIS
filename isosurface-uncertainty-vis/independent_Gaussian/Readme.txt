- Update worklet ContourVariance.cxx to call proper ilerp uncertainty function (Monte Carlo, closed)

- Update the grid resolution in worklet depending on the dataset (we want to automate this in the future)

- Update the isovalue for the contour in run-contour.cxx

- Build viskores freshly using the vtkm source code

- create build directory where you want to build ilerp uncertainty code and "cd" to it

- ccmake source_dir

- For viskores dir point to viskoresBuild/lib/cmake/viskores

- Link with boost, eigen libraries. See CmakeLists.txt in source code

- generate make file and and create executable

- ./run-contour ../correlatedGaussianField.vtk (GaussianField.vtk has two fields mean and variance for uniform volume)

- visualize output.vtk file in ParaView coloerd by ilerp variance