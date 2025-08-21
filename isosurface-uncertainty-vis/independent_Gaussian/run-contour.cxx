#include <viskores/io/VTKDataSetReader.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/cont/Timer.h>
#include <viskores/cont/Initialize.h>

#include "ContourVariance.h"

#include <iostream>

std::string backend = "serial";

void initBackend(viskores::cont::Timer &timer)
{
    // init the vtkh device
    char const *tmp = getenv("UCV_VISKORES_BACKEND");

    if (tmp == nullptr)
    {
        return;
    }
    else
    {
        backend = std::string(tmp);
        std::cout << "Setting the device with UCV_VISKORES_BACKEND=" << backend << "\n";
        std::cout << "This method is antiquated. Consider using the --vtkm-device command line argument." << std::endl;
    }

    // if (rank == 0)
    //{
    std::cout << "vtkm backend is:" << backend << std::endl;
    //}

    if (backend == "serial")
    {
        viskores::cont::RuntimeDeviceTracker &device_tracker = viskores::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(viskores::cont::DeviceAdapterTagSerial());
        timer.Reset(viskores::cont::DeviceAdapterTagSerial());
    }
    else if (backend == "openmp")
    {
        viskores::cont::RuntimeDeviceTracker &device_tracker = viskores::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(viskores::cont::DeviceAdapterTagOpenMP());
        timer.Reset(viskores::cont::DeviceAdapterTagOpenMP());
    }
    else if (backend == "cuda")
    {
        viskores::cont::RuntimeDeviceTracker &device_tracker = viskores::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(viskores::cont::DeviceAdapterTagCuda());
        timer.Reset(viskores::cont::DeviceAdapterTagCuda());
    }
    else
    {
        std::cerr << " unrecognized backend " << backend << std::endl;
    }
    return;
}


int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "USAGE: " << argv[0] << " <filename>.vtk\n";
    return 1;
  }
    
  // Measure time
  auto opts = viskores::cont::InitializeOptions::DefaultAnyDevice;
  viskores::cont::InitializeResult config = viskores::cont::Initialize(argc, argv, opts);

  viskores::cont::Timer timer{ config.Device };
  initBackend(timer);

  viskores::io::VTKDataSetReader reader(argv[1]);
  viskores::cont::DataSet input = reader.ReadDataSet();

  ContourVariance contourFilter;
  contourFilter.SetMeanField("mean");
  contourFilter.SetVarianceField("variance");
  // tangle
  //contourFilter.SetIsoValue(27.6);
    
  // nyc
  contourFilter.SetIsoValue(2);
  contourFilter.SetMergeDuplicatePoints(true);
  contourFilter.SetAddInterpolationEdgeIds(true);
    
  timer.Start();
  viskores::cont::DataSet contour = contourFilter.Execute(input);
  timer.Stop();
  std::cout << "Running time: " << timer.GetElapsedTime() << std::endl;
  

  viskores::io::VTKDataSetWriter writer("output.vtk");
  writer.SetFileTypeToBinary();
  writer.WriteDataSet(contour);

  return 0;
}
