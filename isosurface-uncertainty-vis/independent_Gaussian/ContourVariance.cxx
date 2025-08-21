// Questions:
// Passing rhoX, rhoY, rhoZ
// Passing full ensemble/ array at a voxel to worklet will very useful, if we could do that
// Compiling with openMP. Any changes to cMakeList.txt?
// interpMean is isovalue?

#include "ContourVariance.h"

#include <viskores/cont/ArrayCopy.h>
#include <viskores/Math.h>
#include <viskores/worklet/WorkletMapField.h>

namespace
{

struct ComputeEdgeVarianceWorklet : viskores::worklet::WorkletMapField
{
  using ControlSignature = void(FieldIn edgeIds,
                                FieldIn interpMean,
                                WholeArrayIn inputMeans,
                                WholeArrayIn inputVariance,
                                FieldOut outputVariance);
  using ExecutionSignature = _5(_1, _2, _3, _4);

  template <typename T, typename PortalType>
  VISKORES_EXEC T operator()(const viskores::Id2& edgeIds,
                         T interpMean,
                         const PortalType& inputMeans,
                         const PortalType& inputVariance) const
  {
    // Unfortunately, the contour filter does not currently save the interpolation
    // weight. Thus, we have to estimate it based on the interpolated mean. In the
    // future, it would probably be helpful to optionally save this information
    // because this is silly extra work.
    T weight;
    T mean0;
    T mean1;
    {
      mean0 = inputMeans.Get(edgeIds[0]);
      mean1 = inputMeans.Get(edgeIds[1]);
      T diff = mean1 - mean0;
      if (viskores::Abs(interpMean * 0.000001) > viskores::Abs(diff))
      {
        // Means are basically the same. Avoid divides by zero.
        weight = 0.5;
      }
      else
      {
        weight = (interpMean - mean0) / diff;
      }
    }

    T variance0 = inputVariance.Get(edgeIds[0]);
    T variance1 = inputVariance.Get(edgeIds[1]);

    // Tushar, at this point, I think you should have all the information you
    // need to compute the variance on the edge. I am just doing a simple linear
    // interpolation of the variance, which I know is wrong. I don't understand
    // the math functions you sent, so you'll have to replace this with the
    // correct computation.
      
    gaussian_distribution g;
    viskores::Float64 expt;
    viskores::Float64 cross_prob;
    viskores::Float64 var;
    viskores::Id numSamples = 1000;
    g.gaussian_alpha_pdf(mean0, variance0, mean1, variance1, 0, interpMean, &expt, &cross_prob, &var);
    //g.gaussian_alpha_pdf_MonteCarlo(mean0, variance0, mean1, variance1, 0, interpMean, &expt, &cross_prob, &var, numSamples);
    return var;
    //return vtkm::Lerp(variance0, variance1, weight);
  }
};

template <typename VarianceArrayType>
VISKORES_CONT viskores::cont::UnknownArrayHandle ComputeEdgeVariance(
  const VarianceArrayType& variance,
  const viskores::cont::UnknownArrayHandle& outputMeanArrayUnknown,
  const viskores::cont::UnknownArrayHandle& inputMeanArrayUnknown,
  const viskores::cont::UnknownArrayHandle& edgeIdsUnknown)
{
  VarianceArrayType outputMeanArray;
  viskores::cont::ArrayCopyShallowIfPossible(outputMeanArrayUnknown, outputMeanArray);
  VarianceArrayType inputMeanArray;
  viskores::cont::ArrayCopyShallowIfPossible(inputMeanArrayUnknown, inputMeanArray);
  
  viskores::cont::ArrayHandle<viskores::Id2> edgeIds;
  edgeIdsUnknown.AsArrayHandle(edgeIds);

  viskores::cont::Invoker invoke;
  VarianceArrayType edgeVariance;
  invoke(
    ComputeEdgeVarianceWorklet{}, edgeIds, outputMeanArray, inputMeanArray, variance, edgeVariance);
  return edgeVariance;
}

} // anonymous namespace

viskores::cont::DataSet ContourVariance::DoExecute(const viskores::cont::DataSet& input)
{
  // Do the actual contour extract while saving the edge information.
  viskores::filter::contour::Contour contourExtract;
  contourExtract.SetIsoValues(this->IsoValues);
  contourExtract.SetGenerateNormals(this->GetGenerateNormals());
  contourExtract.SetComputeFastNormals(this->GetComputeFastNormals());
  contourExtract.SetNormalArrayName(this->GetNormalArrayName());
  contourExtract.SetMergeDuplicatePoints(this->GetMergeDuplicatePoints());
  contourExtract.SetActiveField(0, this->GetActiveFieldName(0), this->GetActiveFieldAssociation());

  viskores::filter::FieldSelection fieldSelection = this->GetFieldsToPass();
  // Don't need to interpolate variance information.
  fieldSelection.AddField(this->GetActiveFieldName(1),
                          this->GetActiveFieldAssociation(1),
                          viskores::filter::FieldSelection::Mode::Exclude);
  // Currently need mean information to figure out interpolation weight.
  fieldSelection.AddField(this->GetActiveFieldName(0),
                          this->GetActiveFieldAssociation(0),
                          viskores::filter::FieldSelection::Mode::Select);
  contourExtract.SetFieldsToPass(fieldSelection);

  // Save edge information
  contourExtract.SetAddInterpolationEdgeIds(true);

  viskores::cont::DataSet contours = contourExtract.Execute(input);

  viskores::cont::UnknownArrayHandle outputVariance;
  auto resolveArray = [&](auto varianceArray)
  {
    outputVariance = ComputeEdgeVariance(varianceArray,
                                         this->GetFieldFromDataSet(0, contours).GetData(),
                                         this->GetFieldFromDataSet(0, input).GetData(),
                                         contours.GetField("edgeIds").GetData());
     
  };
  this->CastAndCallScalarField(this->GetFieldFromDataSet(1, input), resolveArray);

  // Add the computed variance.
  contours.AddPointField(this->GetActiveFieldName(1), outputVariance);

  return contours;
}
