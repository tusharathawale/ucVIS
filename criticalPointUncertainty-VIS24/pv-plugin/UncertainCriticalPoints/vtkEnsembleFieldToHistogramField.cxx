#include "vtkEnsembleFieldToHistogramField.h"

#include <EnsembleFieldToHistogramField.h>

#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include <vtkm/cont/Timer.h>

#include "VTKToVTKm.h"
#include "vtkmlib/DataArrayConverters.h"

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkEnsembleFieldToHistogramField);

vtkEnsembleFieldToHistogramField::vtkEnsembleFieldToHistogramField() = default;
vtkEnsembleFieldToHistogramField::~vtkEnsembleFieldToHistogramField() = default;

std::string vtkEnsembleFieldToHistogramField::GetInputArrayName(
  int index, vtkInformationVector** inputVector)
{
  int association = this->GetInputArrayAssociation(index, inputVector);
  vtkDataArray* inputArray = this->GetInputArrayToProcess(index, inputVector);
  if ((association != vtkDataObject::FIELD_ASSOCIATION_POINTS) || (inputArray == nullptr))
  {
    vtkErrorMacro("Invalid array; array missing or not a point array.");
    return 0;
  }

  const char* scalarFieldName = inputArray->GetName();
  if (!scalarFieldName || scalarFieldName[0] == '\0')
  {
    scalarFieldName = tovtkm::NoNameVTKFieldName();
  }

  return scalarFieldName;
}

int vtkEnsembleFieldToHistogramField::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkDataSet* input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet* output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  try
  {
    // Convert the input dataset to VTK-m
    vtkm::cont::DataSet in = vtktovtkm::Convert(input);

    vtkm::filter::uncertainty::EnsembleFieldToHistogramField filter;

    filter.SetNumberOfBins(this->NumberOfBins);

    filter.SetActiveField(this->GetInputArrayName(0, inputVector));

    filter.SetHistogramDensityName(this->HistogramDensityName);
    filter.SetHistogramEdgesName(this->HistogramEdgesName);

    vtkm::cont::Timer timer;
    timer.Start();
    vtkm::cont::DataSet result = filter.Execute(in);
    timer.Stop();
    std::cout << "execution time of ensemble to histogram conversion is " << timer.GetElapsedTime() << " seconds" << std::endl;

    // Convert the result back.
    // It would be easier if there was a simple method to just convert from general
    // vtkm::cont::DataSet to vtkDataSet. However, that does not exist so you have
    // to copy the vtkDataSet structure in VTK and copy the new fields over. I think
    // it is done this way to prevent creating new arrays for what should be shallow
    // copies. (Maybe in the future the data sharing will be good enough where that
    // is not an issue.)
    output->ShallowCopy(input);
    auto copyField = [&](const std::string& fieldName)
    {
      //vtkDataArray* resultingArray = fromvtkm::Convert(result.GetPointField(fieldName));
      vtkDataArray* resultingArray = vtkmtovtk::Convert(result.GetPointField(fieldName));
      if (resultingArray == nullptr)
      {
        vtkWarningMacro(<< "Unable to convert result array " << fieldName << " from VTK-m to VTK.");
        return;
      }
      output->GetPointData()->AddArray(resultingArray);
      resultingArray->FastDelete();
    };
    copyField(this->HistogramDensityName);
    copyField(this->HistogramEdgesName);
  }
  catch (const vtkm::cont::Error& e)
  {
    vtkErrorMacro(<< "VTK-m error: " << e.GetMessage());
    return 0;
  }

  return 1;
}

void vtkEnsembleFieldToHistogramField::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfBins: " << this->NumberOfBins << "\n";
  os << indent << "HistogramDensityName: " << this->HistogramDensityName << "\n";
  os << indent << "HistogramEdgesName: " << this->HistogramEdgesName << "\n";
}

VTK_ABI_NAMESPACE_END
