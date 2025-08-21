#ifndef UCV_ENTROPY_ADAPTIVE_EIGENS_h
#define UCV_ENTROPY_ADAPTIVE_EIGENS_h

#include <viskores/worklet/WorkletMapTopology.h>
#include <cmath>
#include "linalg/EasyLinAlg/eigen.h"

#if defined(VISKORES_CUDA) || defined(VISKORES_KOKKOS_HIP)
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#else
// using the std library
#include <random>
#endif // VISKORES_CUDA

template <int NumVertecies, int NumCases>
class EntropyAdaptiveEigens : public viskores::worklet::WorkletVisitCellsWithPoints
{
public:
    EntropyAdaptiveEigens(double isovalue, int numSamples, double thresholdInd, double thresholdEnergy)
        : m_isovalue(isovalue), m_numSamples(numSamples), m_thresholdInd(thresholdInd), m_thresholdEnergy(thresholdEnergy){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);

    // using ExecutionSignature = void(_2, _3, _4, _5, _6, WorkIndex);
    using ExecutionSignature = void(_2, _3, _4, _5, _6);
    //  the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename InPointFieldVecMean,
              typename OutCellFieldType1,
              typename OutCellFieldType2,
              typename OutCellFieldType3>
    // VTKM_EXEC void operator()(
    //     const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
    //     const InPointFieldVecMean &inMeanArray,
    //     OutCellFieldType1 &outCellFieldCProb,
    //     OutCellFieldType2 &outCellFieldNumNonzeroProb,
    //     OutCellFieldType3 &outCellFieldEntropy, vtkm::Id workIndex) const
    VISKORES_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        const InPointFieldVecMean &inMeanArray,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy) const
    {
        // printf("debug workIndex %lld\n",workIndex);

        // how to process the case where there are multiple variables
        viskores::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();
        // constexpr uint8_t numVertex3d = 8;
        // only support NumVertecies equals to 4 or 8
        if (numVertexies != NumVertecies)
        {
            printf("2d version has 4 vertecies 16(2^4) cases, 3d version has 8 vertecies 256(2^8) cases\n");
            return;
        }

        if (inMeanArray.GetNumberOfComponents() != NumVertecies)
        {
            printf("inMeanArray in 2d version expect 4 vertecies 3d version expects 8 vertecies\n");
            return;
        }

        // set the trim options to filter out values that does not contain the iso value
        // there is no cross prob for this values
        // find min and cell for all cell values
        using VecType = decltype(inPointFieldVecEnsemble[0]);
        double cellMin = viskores::Infinity64();
        double cellMax = viskores::NegativeInfinity64();
        for (int i = 0; i < NumVertecies; i++)
        {
            find_min_max<VecType>(inPointFieldVecEnsemble[i], cellMin, cellMax);
        }

        // printf("---debug workindex %d\n min %lf max %lf\n",workIndex,cellMin,cellMax);

        if (this->m_isovalue < cellMin || this->m_isovalue > cellMax)
        {
            outCellFieldCProb = 0;
            return;
        }

        // std::vector<double> cov_matrix;
        // compute number of values in comatrix
        constexpr viskores::IdComponent covMatrixSize = (NumVertecies + 1) * NumVertecies / 2;
        viskores::Vec<viskores::FloatDefault, covMatrixSize> cov_matrix;
        viskores::IdComponent index = 0;
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                float cov = find_covariance<VecType>(inPointFieldVecEnsemble[p], inPointFieldVecEnsemble[q], inMeanArray[p], inMeanArray[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        // generate sample
        // UCVMATH::vec_t ucvmeanv;
        EASYLINALG::Vec<double, NumVertecies> ucvmeanv;

        for (int i = 0; i < NumVertecies; i++)
        {
            ucvmeanv[i] = inMeanArray[i];
        }

        // generate mean and cov matrix
        // UCVMATH::mat_t ucvcov;
        int covindex = 0;
        EASYLINALG::Matrix<double, NumVertecies, NumVertecies> ucvcov;

        for (int p = 0; p < NumVertecies; ++p)
        {
            for (int q = p; q < NumVertecies; ++q)
            {
                // use the elements at the top half
                // printf("%f ", cov_matrix[covindex]);
                ucvcov[p][q] = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    ucvcov[q][p] = ucvcov[p][q];
                }
                covindex++;
            }
        }

        viskores::IdComponent numSamples = this->m_numSamples;

        // EASYLINALG::Matrix<double, numVertex3d, numVertex3d> A = EASYLINALG::SymmEigenDecomposition(ucvcov8by8, this->m_tolerance, this->m_iterations);
        //  Transform the iso value
        EASYLINALG::Vec<double, NumVertecies> transformIso(0);
        for (int i = 0; i < NumVertecies; i++)
        {
            transformIso[i] = this->m_isovalue - ucvmeanv[i];
        }

        // Compute eigen values
        EASYLINALG::Vec<double, NumVertecies> eigenValues;
        EASYLINALG::SymmEigenValues(ucvcov, this->m_tolerance, this->m_iterations, eigenValues);
        using EigenValuesType = decltype(eigenValues);

        // sorting eigen values
        Sort<EigenValuesType>(eigenValues);

        // make sure eigen values are in descending order
        if (checkOrder<EigenValuesType>(eigenValues) == false)
        {
            // print out
            printf("wrong eigne value sequence\n");
            eigenValues.Show();
        }

        // if the last eigen value is comparatively large
        // all eigen values are important
        // constexpr vtkm::Id lastEigenIndex = NumVertecies-1;
        // if(eigenValues[lastEigenIndex]>this->m_thresholdInd*eigenValues[0]){
        // using indepedent gaussian
        //}else{
        // using mc sampling
        // finding number of important eigen values
        //}

        // Choosing number of important eigen values
        EASYLINALG::Vec<double, NumVertecies> eigenValuesFiltered(0);
        // using first eigen as default case
        int filteredEigenCount = 1;
        eigenValuesFiltered[0] = eigenValues[0];
        // filter out eigen values when it is less then the threshold
        // and not all eigen value are used
        // assuming this->m_thresholdEnergy is larger than 0
        for (int i = 1; i < NumVertecies; i++)
        {
            // TODO, there are 0 eigen values in the dataset
            // even if we set eigengy as small value, we still filter out a lot of them
            if (eigenValues[i] > this->m_thresholdEnergy * eigenValues[0])
            {
                eigenValuesFiltered[filteredEigenCount] = eigenValues[i];
                filteredEigenCount++;
            }
        }

        // if (workIndex == 11282)
        // {
        //     printf("filteredEigenCount is %d\n", filteredEigenCount);
        //     eigenValuesFiltered.Show();
        //     printf("cov matrix is\n");
        //     ucvcov.Show();
        //     printf("ucvmeanv is\n");
        //     ucvmeanv.Show();
        // }

        // printf("debug filteredEigenCount %d\n",filteredEigenCount);

        // Compute eigen vectors only for important eigen values
        EASYLINALG::Vec<EASYLINALG::Vec<double, NumVertecies>, NumVertecies> eigenVectors;
        // how many eigen vector we want to use

        for (int i = 0; i < filteredEigenCount; i++)
        {
            eigenVectors[i] = EASYLINALG::ComputeEigenVectors(ucvcov, eigenValuesFiltered[i], this->m_iterations);
        }

        // if (workIndex == 11282)
        // {
        //     printf("debug eigen vec\n");
        //     for (int i = 0; i < filteredEigenCount; i++)
        //     {
        //         eigenVectors[i].Show();
        //     }
        // }

        viskores::Vec<viskores::FloatDefault, NumCases> probHistogram;

        EASYLINALG::Vec<double, NumVertecies> sample_v;

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm(0, 1);
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm(0, 1);
#endif // VTKM_CUDA

        // init to 0
        for (int i = 0; i < NumCases; i++)
        {
            probHistogram[i] = 0.0;
        }

        for (viskores::Id n = 0; n < numSamples; ++n)
        {
            EASYLINALG::Vec<double, NumVertecies> sampleResults(0);

            for (int i = 0; i < filteredEigenCount; i++)
            {

                // vtkm will return nan for negative sqrt value
                // just filter it out in the previous step when filter
                // the eigen value
                sample_v[i] = viskores::Sqrt(eigenValues[i]) * norm(rng);
            }

            // compute sampled results
            // for each sampled results
            for (int i = 0; i < NumVertecies; i++)
            {
                for (int j = 0; j < filteredEigenCount; j++)
                {
                    // eigen vector of jth eigen value, jth sample element
                    // ith componnet in the eigen vector
                    sampleResults[i] += eigenVectors[j][i] * sample_v[j];
                }
            }

            uint caseValue = 0;
            for (uint i = 0; i < NumVertecies; i++)
            {
                // setting associated position to 1 if iso larger then specific cases
                if (transformIso[i] >= sampleResults[i])
                {
                    caseValue = (1 << i) | caseValue;
                }
            }

            // the associated pos is 0 otherwise
            probHistogram[caseValue] = probHistogram[caseValue] + 1.0;
        }

        // go through probHistogram and compute pro
        for (int i = 0; i < NumCases; i++)
        {
            probHistogram[i] = (probHistogram[i] / (1.0 * numSamples));
            // printf("debug caseValue %d probHistogram %f\n", i, probHistogram[i]);
        }

        // cross probability
        // outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
        outCellFieldCProb = 1.0 - (probHistogram[0] + probHistogram[NumCases - 1]);

        viskores::Id nonzeroCases = 0;
        viskores::FloatDefault entropyValue = 0;
        viskores::FloatDefault templog = 0;
        // compute number of nonzero cases
        // compute entropy
        for (int i = 0; i < NumCases; i++)
        {
            if (probHistogram[i] > 0.0001)
            {
                nonzeroCases++;
                templog = viskores::Log2(probHistogram[i]);
                // if (i != 0 && i != totalNumCases - 1)
                //{
                //     totalnonzeroProb += probHistogram[i];
                // }
            }
            // do not update entropy if the pro is zero
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }

        outCellFieldNumNonzeroProb = nonzeroCases;
        outCellFieldEntropy = entropyValue;

        // if (workIndex == 11282)
        // {
        //     for (int i = 0; i < NumCases; i++)
        //     {
        //         printf("probHis %d is %f\n", i, probHistogram[i]);
        //     }
        //     printf("nonzeroCases is %lld\n", nonzeroCases);
        //     printf("outCellFieldCProb is %f\n", outCellFieldCProb);
        //     printf("outCellFieldEntropy is %f\n", outCellFieldEntropy);
        //     printf("debug worklet index %lld\n", workIndex);
        // }

        // check if eigen value is large to small
        // filter out the eigen values that is less then a specific threshold

        /*
        EASYLINALG::Vec<double, NumVertecies> eigenValuesFiltered(0);
        // use all eigen values when the use_all_eigen is true
        int filteredEigenCount = NumVertecies;

        //TODO, checking the eigen values are in desending sequence
        // decide two cases
        //either adopt the indepedent assumption
        //or the case that decide how many eigen values should be used

        if (this->m_use_all_eigen == true)
        {
            eigenValuesFiltered = eigenValues;
        }
        else
        {
            filteredEigenCount = 0;
            // filter out eigen values when it is less then the threshold
            // and not all eigen value are used
            for (int i = 0; i < NumVertecies; i++)
            {
                if (eigenValues[i] > this->m_eigen_threshold)
                {
                    eigenValuesFiltered[filteredEigenCount] = eigenValues[i];
                    filteredEigenCount++;
                }
            }
        }

        // Compute eigen vectors
        EASYLINALG::Vec<EASYLINALG::Vec<double, NumVertecies>, NumVertecies> eigenVectors;
        // how many eigen vector we want to use

        for (int i = 0; i < filteredEigenCount; i++)
        {
            eigenVectors[i] = EASYLINALG::ComputeEigenVectors(ucvcov, eigenValuesFiltered[i], this->m_iterations);
        }

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm;
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm;
#endif // VTKM_CUDA

        vtkm::Vec<vtkm::FloatDefault, NumCases> probHistogram;

        EASYLINALG::Vec<double, numVertex3d> sample_v;

        // init to 0
        for (int i = 0; i < NumCases; i++)
        {
            probHistogram[i] = 0.0;
        }

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            EASYLINALG::Vec<double, numVertex3d> sampleResults(0);

            for (int i = 0; i < filteredEigenCount; i++)
            {
                // sample_v[i]=vtkm::Sqrt(eigenValues[i])*norm(rng);
                std::normal_distribution<double> norm(0, vtkm::Sqrt(eigenValuesFiltered[i]));
                sample_v[i] = norm(rng);
            }

            // compute sampled results
            // for each sampled results
            for (int i = 0; i < numVertex3d; i++)
            {
                for (int j = 0; j < filteredEigenCount; j++)
                {
                    sampleResults[i] += eigenVectors[j][i] * sample_v[j];
                }
            }

            // go through 8 cases
            uint caseValue = 0;
            for (uint i = 0; i < numVertex3d; i++)
            {
                // setting associated position to 1 if iso larger then specific cases
                if (transformIso[i] >= sampleResults[i])
                {
                    caseValue = (1 << i) | caseValue;
                }
            }

            // the associated pos is 0 otherwise
            probHistogram[caseValue] = probHistogram[caseValue] + 1.0;
        }

        // go through probHistogram and compute pro
        for (int i = 0; i < NumCases; i++)
        {
            probHistogram[i] = (probHistogram[i] / (1.0 * numSamples));
            // printf("debug caseValue %d probHistogram %f\n", i, probHistogram[i]);
        }

        // cross probability
        // outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
        outCellFieldCProb = 1.0 - (probHistogram[0] + probHistogram[255]);

        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault entropyValue = 0;
        vtkm::FloatDefault templog = 0;
        // compute number of nonzero cases
        // compute entropy
        for (int i = 0; i < NumCases; i++)
        {
            if (probHistogram[i] > 0.0001)
            {
                nonzeroCases++;
                templog = vtkm::Log2(probHistogram[i]);
                // if (i != 0 && i != totalNumCases - 1)
                //{
                //     totalnonzeroProb += probHistogram[i];
                // }
            }
            // do not update entropy if the pro is zero
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }

        outCellFieldNumNonzeroProb = nonzeroCases;
        outCellFieldEntropy = entropyValue;
         */
    }

    template <typename VecType>
    VISKORES_EXEC inline void Sort(VecType &arr) const
    {
        viskores::Id num = arr.NUM_COMPONENTS;
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < num - i - 1; j++)
            {
                // compare element i and j
                if (arr[j] < arr[j + 1])
                {
                    // swap
                    auto temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
        return;
    }

    template <typename VecType>
    VISKORES_EXEC bool checkOrder(const VecType &eigenValues) const
    {
        viskores::Id num = eigenValues.NUM_COMPONENTS;
        for (viskores::Id i = 0; i < num - 1; i++)
        {
            if (eigenValues[i] < eigenValues[i + 1])
            {
                // in ascending order
                return false;
            }
        }
        return true;
    }

    template <typename VecType>
    VISKORES_EXEC void find_min_max(const VecType &arr, viskores::Float64 &min, viskores::Float64 &max) const
    {
        viskores::Id num = arr.GetNumberOfComponents();
        for (viskores::Id i = 0; i < num; i++)
        {
            // the second one is runtime thing (recombine vec), so convert it into float defualt firt
            viskores::Float64 v = arr[i];
            min = viskores::Min(min, v);
            max = viskores::Max(max, v);
        }
        return;
    }

    template <typename VecType>
    VISKORES_EXEC inline viskores::FloatDefault find_covariance(const VecType &arr1, const VecType &arr2,
                                                        const viskores::FloatDefault &mean1, const viskores::FloatDefault &mean2) const
    {
        if (arr1.GetNumberOfComponents() != arr2.GetNumberOfComponents())
        {
            printf("error, failed to compute find_covariance, the array size should be equal with each other\n");
            return 0;
        }
        viskores::Id arraySize = arr1.GetNumberOfComponents();
        viskores::FloatDefault sum = 0;
        for (int i = 0; i < arraySize; i++)
        {
            viskores::FloatDefault v1 = arr1[i];
            viskores::FloatDefault v2 = arr2[i];
            sum = sum + (v1 - mean1) * (v2 - mean2);
        }

        return sum / (viskores::FloatDefault)(arraySize - 1);
    }

private:
    double m_isovalue;
    int m_numSamples;
    int m_iterations = 200;
    double m_tolerance = 0.0001;

    // threshold to depend if there is sphere covaraince structure
    double m_thresholdInd = 1.0;
    // threshold to determine the number of eigen values we should keep
    double m_thresholdEnergy = 0.1;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN3D_h
