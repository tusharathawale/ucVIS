#ifndef UCV_EXTRACTING_MEAN_STD_h
#define UCV_EXTRACTING_MEAN_STD_h

#include <viskores/worklet/WorkletMapField.h>
#include <viskores/worklet/WorkletReduceByKey.h>
#include <cmath>

struct ExtractingMean : public viskores::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldOut);
    using ExecutionSignature = void(_1, _2);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VISKORES_EXEC void operator()(
        const OriginalValuesType &inPointFieldVecEnsemble, OutputType &meanValue) const
    {
        viskores::FloatDefault boxSum = 0;

        viskores::IdComponent NumComponents = inPointFieldVecEnsemble.GetNumberOfComponents();

        // refer to https://www.strchr.com/standard_deviation_in_one_pass
        for (viskores::IdComponent index = 0;
             index < NumComponents; index++)
        {
            boxSum = boxSum + static_cast<viskores::FloatDefault>(inPointFieldVecEnsemble[index]);
        }

        meanValue = boxSum / (1.0 * (NumComponents));
    }
};

// go through each vertexies and compute associated mean and stdev
struct ExtractingMeanStdevEnsembles : public viskores::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VISKORES_EXEC void operator()(
        const OriginalValuesType &inPointFieldVecEnsemble, OutputType &meanValue, OutputType &stdevValue) const
    {
        viskores::FloatDefault boxSum = 0;

        viskores::IdComponent NumComponents = inPointFieldVecEnsemble.GetNumberOfComponents();

        // refer to https://www.strchr.com/standard_deviation_in_one_pass
        for (viskores::IdComponent index = 0;
             index < NumComponents; index++)
        {
            boxSum = boxSum + static_cast<viskores::FloatDefault>(inPointFieldVecEnsemble[index]);
        }

        meanValue = boxSum / (1.0 * (NumComponents));

        viskores::FloatDefault diffSum = 0;

        for (viskores::IdComponent index = 0;
             index < NumComponents; index++)
        {
            viskores::FloatDefault diff = static_cast<viskores::FloatDefault>(inPointFieldVecEnsemble[index]) - static_cast<viskores::FloatDefault>(meanValue);
            diffSum += diff * diff;
        }

        stdevValue = std::sqrt(diffSum / (1.0*NumComponents));
    }
};

struct ExtractingMeanStdev : public viskores::worklet::WorkletReduceByKey
{
    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut, ReducedValuesOut);
    using ExecutionSignature = void(_2, _3, _4);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VISKORES_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &meanValue, OutputType &stdevValue) const
    {
        viskores::FloatDefault boxSum = 0;

        viskores::IdComponent NumComponents = originalValues.GetNumberOfComponents();

        // refer to https://www.strchr.com/standard_deviation_in_one_pass
        for (viskores::IdComponent index = 0;
             index < NumComponents; index++)
        {
            boxSum = boxSum + static_cast<viskores::FloatDefault>(originalValues[index]);
        }

        meanValue = boxSum / (1.0 * (NumComponents));

        viskores::FloatDefault diffSum = 0;

        for (viskores::IdComponent index = 0;
             index < NumComponents; index++)
        {
            viskores::FloatDefault diff = static_cast<viskores::FloatDefault>(originalValues[index]) - static_cast<viskores::FloatDefault>(meanValue);
            diffSum += diff * diff;
        }

        stdevValue = std::sqrt(diffSum / (1.0*NumComponents));
    }
};

#endif // UCV_EXTRACTING_MEAD_STD_h
