#include <viskores/filter/contour/Contour.h>
#include "gaussian_distribution.h"
#include "uniform_even_z_density.h"

class ContourVariance : public viskores::filter::contour::AbstractContour
{
public:
  VISKORES_CONT void SetMeanField(const std::string& name)
  {
    this->SetActiveField(0, name, viskores::cont::Field::Association::Points);
  }
  VISKORES_CONT void SetVarianceField(const std::string& name)
  {
    this->SetActiveField(1, name, viskores::cont::Field::Association::Points);
  }

protected:
  VISKORES_CONT viskores::cont::DataSet DoExecute(const viskores::cont::DataSet &input) override;
};
