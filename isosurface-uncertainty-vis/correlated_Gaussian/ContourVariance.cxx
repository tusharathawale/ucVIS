// Questions:
// Passing rhoX, rhoY, rhoZ
// Passing full ensemble/ array at a voxel to worklet will very useful, if we could do that
// Compiling with openMP. Any changes to cMakeList.txt?
// interpMean is isovalue?

#include "ContourVariance.h"

#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/CellSetStructured.h>
#include <viskores/Math.h>
#include <viskores/worklet/WorkletMapField.h>

namespace
{

struct ComputeEdgeVarianceWorklet : viskores::worklet::WorkletMapField
{
  using ControlSignature = void(FieldIn edgeIds,
                                FieldIn interpMean,
                                WholeArrayIn inputMeans,
                                WholeArrayIn inputVariance, WholeArrayIn inputRhoX, WholeArrayIn inputRhoY, WholeArrayIn inputRhoZ,
                                FieldOut outputVariance);
  using ExecutionSignature = _8(_1, _2, _3, _4, _5, _6, _7);

  VISKORES_CONT ComputeEdgeVarianceWorklet(const viskores::Id3 resolution)
    : Resolution(resolution)
  {}

  template <typename T, typename PortalType>
  VISKORES_EXEC T operator()(const viskores::Id2& edgeIds,
                         T interpMean,
                         const PortalType& inputMeans,
                         const PortalType& inputVariance, const PortalType& inputRhoX, const PortalType& inputRhoY, const PortalType& inputRhoZ) const
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
    // Manually set grid resolution per dataset
        
    // Convert edgeID to indX, indY, indZ based of grid resolution
      
    /*viskores::Id a0 = edgeIds[0]/(this->Resolution[1]*this->Resolution[2]);
    viskores::Id b0 = (edgeIds[0]%(this->Resolution[1]*this->Resolution[2]))/this->Resolution[2];
    viskores::Id c0 = edgeIds[0]%this->Resolution[2];
          
    viskores::Id a1 = edgeIds[1]/(this->Resolution[1]*this->Resolution[2]);
    viskores::Id b1 = (edgeIds[1]%(this->Resolution[1]*this->Resolution[2]))/this->Resolution[2];
    viskores::Id c1 = edgeIds[1]%this->Resolution[2];*/
      
    // This is correct! and gives smooth uncertainty result for NYC dataset. Might have to play around for other datasets depending on
    // which dimension is referred to 0 and 1 in other datasets
    viskores::Id a0 = edgeIds[0]/(this->Resolution[1]*this->Resolution[0]);
    viskores::Id b0 = (edgeIds[0]%(this->Resolution[1]*this->Resolution[0]))/this->Resolution[0];
    viskores::Id c0 = edgeIds[0]%this->Resolution[0];
            
    viskores::Id a1 = edgeIds[1]/(this->Resolution[1]*this->Resolution[0]);
    viskores::Id b1 = (edgeIds[1]%(this->Resolution[1]*this->Resolution[0]))/this->Resolution[0];
    viskores::Id c1 = edgeIds[1]%this->Resolution[0];

    viskores::Float64 covar=0.0;
      
    // Extract covariance depending on if edge is going in, X, Y, Z direction
    
    //(a,X;b,Y;c,Z)
    if (viskores::Abs(a1 - a0) == 1)
    {
        if(a0 < a1)
            covar = inputRhoX.Get(edgeIds[0]);
        else
            covar = inputRhoX.Get(edgeIds[1]);
    }
        
    else if (viskores::Abs(b1 - b0) == 1)
    {
        if(b0 < b1)
            covar = inputRhoY.Get(edgeIds[0]);
        else
            covar = inputRhoY.Get(edgeIds[1]);
    }
        
    else if (viskores::Abs(c1 - c0) == 1)
    {
        if(c0 < c1)
            covar = inputRhoZ.Get(edgeIds[0]);
        else
            covar = inputRhoZ.Get(edgeIds[1]);
    }
      
   //(a,X;b,Y;c,Z), (a,X;b,Z;c,Y), (a,Y;b,X;c,Z), (a,Y;b,Z;c,X), (a,Z;b,X;c,Y), (a,Z;b,Y;c,X)
    
    // Covariance is multiplication of correlation with standard dev of individual vertices
    //covar = covar*viskores::Sqrt(variance0)*viskores::Sqrt(variance1);
      
    // Currently, rhoX, rhoY, and rhoZ fields store covariance and not correlation, based on the python code (this should be updated to avoid confusion)
    //std::cout<<covar/(viskores::Sqrt(variance0)*viskores::Sqrt(variance1))<<"\n";
      
    // Tushar, at this point, I think you should have all the information you
    // need to compute the variance on the edge. I am just doing a simple linear
    // interpolation of the variance, which I know is wrong. I don't understand
    // the math functions you sent, so you'll have to replace this with the
    // correct computation.
      
    gaussian_distribution g;
    viskores::Float64 expt;
    viskores::Float64 cross_prob;
    viskores::Float64 var;
    viskores::Float64 varTemp;
    viskores::Id numSamples = 4000;
      
    // Closed-form with correlated assumption
    g.gaussian_alpha_pdf(mean0, variance0, mean1, variance1, covar, interpMean, &expt, &cross_prob, &var);
    
    // Monte Carlo with correlated assumption
    //g.gaussian_alpha_pdf_MonteCarlo(mean0, variance0, mean1, variance1, covar, interpMean, &expt, &cross_prob, &var, numSamples);
    //cout<<"var1:"<<variance0<<", var2:"<<variance1<<", covar:"<<covar<<"ilerp var<<"<<var<<"\n";
      
    // Code to calculate absolute difference in variances of independent and correlated Gaussian assumption
    // Closed-form with independent assumption
    //g.gaussian_alpha_pdf(mean0, variance0, mean1, variance1, 0, interpMean, &expt, &cross_prob, &var);
    //g.gaussian_alpha_pdf(mean0, variance0, mean1, variance1, 0, interpMean, &expt, &cross_prob, &varTemp);
      
    //Calculate difference: Comment out the part after use
    //var = viskores::Abs(var-varTemp);
    //var = var-varTemp;
      
    //if (var < 0)
    //    std::cout<<var<<"\n";
      
    return var;
    //return vtkm::Lerp(variance0, variance1, weight);
  }

private:
  viskores::Id3 Resolution;
};

template <typename VarianceArrayType>
VISKORES_CONT viskores::cont::UnknownArrayHandle ComputeEdgeVariance(
  const VarianceArrayType& variance, const viskores::cont::UnknownArrayHandle& inputRhoXArrayUnknown, const viskores::cont::UnknownArrayHandle& inputRhoYArrayUnknown, const viskores::cont::UnknownArrayHandle& inputRhoZArrayUnknown,
  const viskores::cont::UnknownArrayHandle& outputMeanArrayUnknown,
  const viskores::cont::UnknownArrayHandle& inputMeanArrayUnknown,
  const viskores::cont::UnknownArrayHandle& edgeIdsUnknown,
  viskores::Id3 resolution)
{
  VarianceArrayType outputMeanArray;
  viskores::cont::ArrayCopyShallowIfPossible(outputMeanArrayUnknown, outputMeanArray);
    
  VarianceArrayType inputMeanArray;
  viskores::cont::ArrayCopyShallowIfPossible(inputMeanArrayUnknown, inputMeanArray);
    
  VarianceArrayType inputRhoXArray;
  viskores::cont::ArrayCopyShallowIfPossible(inputRhoXArrayUnknown, inputRhoXArray);
  
  VarianceArrayType inputRhoYArray;
  viskores::cont::ArrayCopyShallowIfPossible(inputRhoYArrayUnknown, inputRhoYArray);
    
  VarianceArrayType inputRhoZArray;
  viskores::cont::ArrayCopyShallowIfPossible(inputRhoZArrayUnknown, inputRhoZArray);
  
  viskores::cont::ArrayHandle<viskores::Id2> edgeIds;
  edgeIdsUnknown.AsArrayHandle(edgeIds);
    
  //std::cout<<inputRhoXArray[0];

  viskores::cont::Invoker invoke;
  VarianceArrayType edgeVariance;
  invoke(
    ComputeEdgeVarianceWorklet{ resolution }, edgeIds, outputMeanArray, inputMeanArray, variance, inputRhoXArray, inputRhoYArray, inputRhoZArray, edgeVariance);
  return edgeVariance;
}

} // anonymous namespace

viskores::cont::DataSet ContourVariance::DoExecute(const viskores::cont::DataSet& input)
{
  // Get the resolution of the input (which must be structured) so that we can
  // find the edge data.
  viskores::cont::CellSetStructured<3> cellSet;
  input.GetCellSet().AsCellSet(cellSet);
  viskores::Id3 resolution = cellSet.GetPointDimensions();

  // Do the actual contour extract while saving the edge information.
  viskores::filter::contour::Contour contourExtract;
  contourExtract.SetIsoValues(this->IsoValues);
  contourExtract.SetGenerateNormals(this->GetGenerateNormals());
  contourExtract.SetComputeFastNormals(this->GetComputeFastNormals());
  contourExtract.SetNormalArrayName(this->GetNormalArrayName());
  contourExtract.SetMergeDuplicatePoints(this->GetMergeDuplicatePoints());
  contourExtract.SetActiveField(0, this->GetActiveFieldName(0), this->GetActiveFieldAssociation());

  viskores::filter::FieldSelection fieldSelection = this->GetFieldsToPass();
    
  // Currently need mean information to figure out interpolation weight.
  fieldSelection.AddField(this->GetActiveFieldName(0),
                          this->GetActiveFieldAssociation(0),
                          viskores::filter::FieldSelection::Mode::Select);
    
  // Don't need to interpolate variance information.
  fieldSelection.AddField(this->GetActiveFieldName(1),
                          this->GetActiveFieldAssociation(1),
                          viskores::filter::FieldSelection::Mode::Exclude);
  
  // Don't need to interpolate rhoX information.
  fieldSelection.AddField(this->GetActiveFieldName(2),
                          this->GetActiveFieldAssociation(2),
                          viskores::filter::FieldSelection::Mode::Exclude);
    
  // Don't need to interpolate rhoY information.
  fieldSelection.AddField(this->GetActiveFieldName(3),
                          this->GetActiveFieldAssociation(3),
                          viskores::filter::FieldSelection::Mode::Exclude);
    
  // Don't need to interpolate rhoZ information.
  fieldSelection.AddField(this->GetActiveFieldName(4),
                          this->GetActiveFieldAssociation(4),
                          viskores::filter::FieldSelection::Mode::Exclude);
  
  contourExtract.SetFieldsToPass(fieldSelection);

  // Save edge information
  contourExtract.SetAddInterpolationEdgeIds(true);

  viskores::cont::DataSet contours = contourExtract.Execute(input);

  viskores::cont::UnknownArrayHandle outputVariance;
  auto resolveArray = [&](auto varianceArray)
  {
    // pass variance array, rhoX, rhoY, and rhoZ, mean from contour, mean form input, edge Ids
    outputVariance = ComputeEdgeVariance(varianceArray,
                                         this->GetFieldFromDataSet(2, input).GetData(),
                                         this->GetFieldFromDataSet(3, input).GetData(),
                                         this->GetFieldFromDataSet(4, input).GetData(),
                                         this->GetFieldFromDataSet(0, contours).GetData(),
                                         this->GetFieldFromDataSet(0, input).GetData(),
                                         contours.GetField("edgeIds").GetData(),
                                         resolution);
     
  };
  this->CastAndCallScalarField(this->GetFieldFromDataSet(1, input), resolveArray);

  // Add the computed variance.
  contours.AddPointField(this->GetActiveFieldName(1), outputVariance);

  return contours;
}
