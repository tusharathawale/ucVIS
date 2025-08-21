#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/io/VTKDataSetReader.h>

#include <viskores/cont/Initialize.h>
#include <viskores/cont/ArrayHandleRuntimeVec.h>
#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/ArrayHandleRuntimeVec.h>
#include <viskores/cont/Timer.h>
#include <viskores/cont/DataSetBuilderUniform.h>

// #include "ucvworklet/ExtractMinMaxOfPoint.hpp"
// #include "ucvworklet/CriticalPointMonteCarlo.hpp"
#include <viskores/worklet/WorkletMapField.h>
#include <math.h>
#include <float.h>
#include "MVGaussianTemp.hpp"
#include "EntropyAdaptiveEigens.hpp"
#include "ExtractingMeanStdev.hpp"

using SupportedTypesVec = viskores::List<viskores::Vec<float, 7>>;


viskores::cont::ArrayHandle<viskores::Float64> ComputeEntropyWithRuntimeVec(viskores::cont::DataSet vtkmDataSet,
                                                                    double isovalue, int numSamples, std::string outputFileNameSuffix, bool use2d, double eigenThreshold, viskores::cont::Timer &timer, std::string writeFile)
{
    timer.Start();
    // Processing current ensemble data sets based on uncertianty countour
    viskores::cont::ArrayHandle<viskores::FloatDefault> meanArray;
    viskores::cont::ArrayHandle<viskores::Float64> crossProbability;
    viskores::cont::ArrayHandle<viskores::Id> numNonZeroProb;
    viskores::cont::ArrayHandle<viskores::Float64> entropy;
    auto resolveType = [&](auto &concreteArray)
    {
        viskores::cont::Invoker invoke;
        viskores::Id numPoints = concreteArray.GetNumberOfValues();
        auto concreteArrayView = viskores::cont::make_ArrayHandleView(concreteArray, 0, numPoints);

        invoke(ExtractingMean{}, concreteArrayView, meanArray);
        // printSummary_ArrayHandle(meanArray, std::cout);
        // printSummary_ArrayHandle(stdevArray, std::cout);

        // check 2d or 3d
        
        // check 2d or 3d
        if (use2d)
        {
            invoke(EntropyAdaptiveEigens<4, 16>{isovalue, numSamples, 1.0, eigenThreshold}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }
        else
        {
            invoke(EntropyAdaptiveEigens<8, 256>{isovalue, numSamples, 1.0, eigenThreshold}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }

        auto outputDataSet = vtkmDataSet;
        
        timer.Stop();
        std::cout << "worklet time is:" << timer.GetElapsedTime() << std::endl;

        std::stringstream stream;
        stream << isovalue;
        std::string isostr = stream.str();

        std::string outputFileName = outputFileNameSuffix + isostr + ".vtk";

        outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
        outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
        outputDataSet.AddCellField("entropy" + isostr, entropy);
        if (writeFile == "true")
        {
            viskores::io::VTKDataSetWriter writeCross(outputFileName);
            writeCross.WriteDataSet(outputDataSet);
        }
    };

    vtkmDataSet.GetField("ensembles")
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);

    return crossProbability;
}



void callWorklet(viskores::cont::DataSet& vtkmDataSet, double isovalue, viskores::Id numSamples)
{
    viskores::cont::ArrayHandle<viskores::Float64> crossProbability;
    viskores::cont::ArrayHandle<viskores::Id> numNonZeroProb;
    viskores::cont::ArrayHandle<viskores::Float64> entropy;
 // executing the uncertianty thing
 using WorkletType = MVGaussianTemp;
 using DispatcherType = viskores::worklet::DispatcherMapTopology<WorkletType>;
 auto resolveType = [&](const auto &concrete)
 {
   DispatcherType dispatcher(MVGaussianTemp{isovalue, numSamples});
   dispatcher.Invoke(vtkmDataSet.GetCellSet(), concrete, crossProbability, numNonZeroProb, entropy);
 };
    
std::cout << "Came here! " <<std::endl;
vtkmDataSet.PrintSummary(std::cout);
vtkmDataSet.GetField("ensembles").GetData().CastAndCallForTypes<SupportedTypesVec, VISKORES_DEFAULT_STORAGE_LIST>(resolveType);

std::cout << "Came here after! " <<std::endl;
    
std::stringstream stream;
//stream << std::fixed << std::setprecision(2) << isovalue;
stream << isovalue;
std::string isostr = stream.str();

// check results
// we use a shallow copy as the data set for
auto outputDataSet = vtkmDataSet;
outputDataSet.AddCellField("cross_prob", crossProbability);
outputDataSet.AddCellField("num_nonzero_prob", numNonZeroProb);
outputDataSet.AddCellField("entropy", entropy);

std::string outputFileName = "./ucv_3d_iso_" + isostr + ".vtk";
viskores::io::VTKDataSetWriter writeCross(outputFileName);
writeCross.WriteDataSet(outputDataSet);
    
}


int main(int argc, char *argv[])
{
    viskores::cont::InitializeResult initResult = viskores::cont::Initialize(
        argc, argv, viskores::cont::InitializeOptions::DefaultAnyDevice);
    viskores::cont::Timer timer{initResult.Device};

    if (argc != 12)
    {
        //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
        std::cout << "<executable> <SyntheticDataSuffix> <FieldName> <Dimx> <Dimy> <Dimz> <num of ensembles> <num of samples> <isovalue> <outputFileSuffix> <eigenThreshold> <true/false (write file or not)>" << std::endl;
        exit(0);
    }

    std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

    std::string dataPathSuffix = std::string(argv[1]);
    std::string fieldName = std::string(argv[2]);

    int dimx = std::stoi(argv[3]);
    int dimy = std::stoi(argv[4]);
    int dimz = std::stoi(argv[5]);

    int numEnsembles = std::stoi(argv[6]);
    int numSamples = std::stoi(argv[7]);
    
    double isovalue = std::stod(argv[8]);
    
    std::string outputSuffix = std::string(argv[9]);
    double eigenThreshold = std::stod(argv[10]);
    std::string writeFile = std::string(argv[11]);
    
    if (eigenThreshold < 0)
    {
        throw std::runtime_error("eigenThreshold is supposed to be larger than 0");
    }


    //update spacing and origin as needed by the dataset4
    const viskores::Id3 dims(dimx, dimy, dimz);
    const viskores::Id3 origin(0, 0, 0);
    const viskores::Id3 spacing(2, 2, 2);
    viskores::cont::DataSetBuilderUniform dataSetBuilder;
    viskores::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims,origin,spacing);

    viskores::cont::ArrayHandleRuntimeVec<viskores::FloatDefault> allEnsemblesArray(numEnsembles);
    allEnsemblesArray.Allocate(dimx * dimy * dimz);

    std::vector<viskores::cont::ArrayHandle<viskores::FloatDefault>> dataArray;
    // redsea data start from 1
    for (int ensId = 0; ensId < numEnsembles; ensId++)
    {
        std::string fileName = dataPathSuffix + std::to_string(ensId) + ".vtk";
        viskores::io::VTKDataSetReader reader(fileName);
        viskores::cont::DataSet inData = reader.ReadDataSet();

        viskores::cont::ArrayHandle<viskores::FloatDefault> fieldDataArray;
        viskores::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldName).GetData(), fieldDataArray);
        dataArray.push_back(fieldDataArray);
        //printSummary_ArrayHandle(fieldDataArray, std::cout, true);

    }

    std::cout << "ok to load the data at the first step" << std::endl;

    // using all ensembles
    viskores::cont::ArrayHandleRuntimeVec<viskores::FloatDefault> runtimeVecArray(numEnsembles);
    runtimeVecArray.Allocate(dimx * dimy * dimz);
    auto writePortal = runtimeVecArray.WritePortal();
    for (int k = 0; k < dimz; k++)
    {
        for (int j = 0; j < dimy; j++)
        {
            for (int i = 0; i < dimx; i++)
            {
                int pointIndex = k * dimx * dimy + j * dimx + i;
                //printf("debug pointIndex %d\n",pointIndex);
                auto vecValue = writePortal.Get(pointIndex);
                // load ensemble data from ens 0 to ens with id usedEnsembles-1
                for (int currEndId = 0; currEndId < numEnsembles; currEndId++)
                {
                    // set ensemble value
                    vecValue[currEndId] = dataArray[currEndId].ReadPortal().Get(pointIndex);
                }
            }
        }
    }

    vtkmDataSet.AddPointField("ensembles", runtimeVecArray);
    printSummary_ArrayHandle(runtimeVecArray, std::cout);

    // using pointNeighborhood worklet to process the data
    // timer.Start();
    // callWorklet(vtkmDataSet, isovalue, numSamples);
    // timer.Stop();
    
    
    bool use2d = false;
    if (dimz == 1)
    {
        use2d = true;
    }

    std::cout << "start to call the adaptive worklet" << std::endl;
    auto crossProb1 = ComputeEntropyWithRuntimeVec(vtkmDataSet, isovalue, numSamples, outputSuffix, use2d, eigenThreshold, timer, writeFile);
    
    // std::cout << "filter execution time: " << timer.GetElapsedTime() << std::endl;

    // std::string outputFileName = "MinProb_MC" + std::to_string(dimx) + "_" + std::to_string(dimy) + "ens_" + std::to_string(numEnsembles) + ".vtk";
    // vtkm::io::VTKDataSetWriter writeCross(outputFileName);
    // writeCross.WriteDataSet(vtkmDataSet);

    return 0;
}
