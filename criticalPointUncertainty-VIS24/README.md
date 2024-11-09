This code uses the vtk-m to implement the uncertainty algorithm to find the critical points. The python scripts are provided too in pythonScripts folder.

The details of the critical point uncertainty methods can be found in the following paper:

Uncertainty Visualization of Critical Points of 2D Scalar Fields for Parametric and Nonparametric Probabilistic Models
Tushar M. Athawale, Zhe Wang, Kenneth Moreland, David Pugmire, Qian Gong, Scott Klasky, Chris R. Johnson, Paul Rosen, "Uncertainty Visualization of Critical Points of 2D Scalar Fields for Parametric and Nonparametric Probabilistic Models," in IEEE Transactions on Visualization and Computer Graphics, doi: 10.1109/TVCG.2024.3456393.

### Build testing files

The installing scripts are listed under the `UCV/install_scripts`.

For example, on the mac platform, we need to execute 
`sh mac.sh`. To build plugin, mark plugin flag ON during ccmake and  run `sh mac_plugin_paraview.sh`. Update script files to appropriately point to the path of vtkm and QT5 and vtkm version for manual vtkm installation.

under the install_scripts folder. This script will create a folder called mac. Under this folder, the dependency, namely the vtk-m, will be downloaded under the `src` folder and installed into a specific folder. The executable file will be installed at `UCV/install_scripts/mac/install/UCV`

The `exp_scripts/frontier_run_weather.sh` is the script to run experiments on the Frontier supercomputer.

### Build ParaView plugin

We first need to build the ParaView on a specific platform. The user can refer to [this document](https://gitlab.kitware.com/paraview/paraview/blob/master/Documentation/dev/build.md) to build the ParaView manually.

When installing the uncertainty critical point filter through ParaView, we need to set

`DParaView_DIR` and `DVTKm_DIR` as the proper path (see the provided script mac_plugin_paraview.sh). In this case, we use the vtk-m embedded with ParaView download as the third-party dependency of during ParaView installation.

After installing the ParaView plugin, the user could load the associated plugin through the ParaView filter and run the loaded filter. The plugin is a .so file located in UCV/install_scripts/mac/install/UCV/lib/UncertainCriticalPoints

The [demo video](https://drive.google.com/file/d/1GS0OJW_HQWHP5HyS8xV0cxbDHKK_sRgR/view?usp=sharing) shows how to run the filter through the ParaView plugin.

The input data file to run a plugin comprises a single vtk file that has multiple 2D ensemble members

### Additional Comments
-Build paraview against QT 5. By default QT6 will be configured in ccmake. On our system, Qt5 is installed in /usr/local/Cellar/qt@5/5.15.13_1/lib/cmake/Qt5. This will also ensure that there are no errors when generating make file.

- You may need to run cmake -DPARAVIEW_QT_VERSION=5 -DVTK_QT_VERSION=5 ../src/ and then run ccmake for additional configuration e.g., tbb etc. to resolve conflict with Qt6 and successful makefile generation

-Also be consistent in QT 5 paths. In a generated CMakeCache.txt if Qt 5 are referring to two differently installed folders (e.g., /usr/local/Cellar/qt@5 and /usr/local/anaconda) it will not build paraview successfully. Update the paths to only point to /usr/local/Cellar/qt@5/

-Lastly, during ParaView installation, make sure VTK-m dir is pointed to the one that comes with ParaView Source code and not VTK-m manually downloaded.

