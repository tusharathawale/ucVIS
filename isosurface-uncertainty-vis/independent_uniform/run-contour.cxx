#include <viskores/io/VTKDataSetReader.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/cont/Timer.h>
#include <viskores/cont/Initialize.h>

#include "ContourVariance.h"

#include <iostream>

int main(int argc, char** argv)
{
  // Measure time
  //auto opts = viskores::cont::InitializeOptions::RequireDevice;
  auto opts = viskores::cont::InitializeOptions::DefaultAnyDevice;
  viskores::cont::InitializeResult config = viskores::cont::Initialize(argc, argv, opts);

  if (argc != 2) {
    std::cerr << "USAGE: " << argv[0] << " [Viskores options] <filename>.vtk\n\n";
    std::cerr << "Viskores options are:\n" << config.Usage;
    return 1;
  }

  if (getenv("UCV_VISKORES_BACKEND")) {
    std::cerr << "WARNING: UCV_VISKORES_BACKEND is no longer supported and will be ignored.\n";
  }
  std::cout << "Using Vikores device " << config.Device.GetName() << "\n";

  viskores::cont::Timer timer{ config.Device };

  viskores::io::VTKDataSetReader reader(argv[1]);
  viskores::cont::DataSet input = reader.ReadDataSet();

  ContourVariance contourFilter;
  contourFilter.SetMeanField("mean");
  // For independent uniform distribution, variance is half of the uniform distribution width
  contourFilter.SetVarianceField("variance");
  //contourFilter.SetRhoXField("rhoX");
  //contourFilter.SetRhoYField("rhoY");
  //contourFilter.SetRhoZField("rhoZ");
    
  //tangle
  //contourFilter.SetIsoValue(27.6);
    
  //nyc
  contourFilter.SetIsoValue(2);

  contourFilter.SetMergeDuplicatePoints(true);
  contourFilter.SetAddInterpolationEdgeIds(true);

  // Do "burn-in" run to let device warm up.
  viskores::cont::DataSet contour = contourFilter.Execute(input);

  constexpr viskores::IdComponent NUM_TRIALS = 1;
  viskores::Float64 totalTime = 0;
  for (viskores::IdComponent trial = 0; trial < NUM_TRIALS; ++trial)
  {
    timer.Start();
    contour = contourFilter.Execute(input);
    timer.Stop();
    viskores::Float64 trialTime = timer.GetElapsedTime();
    std::cout << "Running time (seconds): " << trialTime << std::endl;
    totalTime += trialTime;
  }

  std::cout << "\nAverage running time (seconds): " << totalTime/NUM_TRIALS << "\n";

  viskores::io::VTKDataSetWriter writer("output.vtk");
  writer.SetFileTypeToBinary();
  writer.WriteDataSet(contour);

  return 0;
}
