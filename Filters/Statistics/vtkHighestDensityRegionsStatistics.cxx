#include "vtkHighestDensityRegionsStatistics.h"

#include "vtkDataArrayCollection.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkStatisticsAlgorithmPrivate.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"

#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <sstream>

vtkStandardNewMacro(vtkHighestDensityRegionsStatistics);

// ----------------------------------------------------------------------
vtkHighestDensityRegionsStatistics::vtkHighestDensityRegionsStatistics()
  {
  // Initialize H smooth matrix to Identity.
  this->SmoothHC1[0] = 1.0;
  this->SmoothHC1[1] = 0.0;
  this->SmoothHC2[0] = 0.0;
  this->SmoothHC2[1] = 1.0;

  //  At the construction, no columns pair are requested yet
  this->NumberOfRequestedColumnsPair = 0;
  }

// ----------------------------------------------------------------------
vtkHighestDensityRegionsStatistics::~vtkHighestDensityRegionsStatistics()
  {
  }

// ----------------------------------------------------------------------
void vtkHighestDensityRegionsStatistics::PrintSelf(ostream& os,
                                                   vtkIndent indent)
  {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Smooth matrix: " <<
    this->SmoothHC1[0] << ", " <<
    this->SmoothHC1[1] << ", " <<
    this->SmoothHC2[0] << ", " <<
    this->SmoothHC2[1] << "\n";
  }

// ----------------------------------------------------------------------
void vtkHighestDensityRegionsStatistics::SetSigma(double sigma)
  {
  if (this->SmoothHC1[0] == sigma &&
    this->SmoothHC1[1] == 0.0 &&
    this->SmoothHC2[0] == 0.0 &&
    this->SmoothHC2[1] == sigma)
    {
    return;
    }
  // Force H matrix to be equal to sigma * Identity.
  this->SmoothHC1[0] = sigma;
  this->SmoothHC1[1] = 0.0;
  this->SmoothHC2[0] = 0.0;
  this->SmoothHC2[1] = sigma;
  this->Modified();
  }

// ----------------------------------------------------------------------
void vtkHighestDensityRegionsStatistics::Learn(vtkTable* inData,
                                               vtkTable* vtkNotUsed(inParameters),
                                               vtkMultiBlockDataSet* outMeta)
  {
  if (!inData || !outMeta)
    {
    return;
    }

  vtkNew<vtkTable> outputColumns;

  std::set<std::set<vtkStdString> >::const_iterator reqIt;

  // Make sure the number of requested pairs of columns is 0
  // before the computation.
  this->NumberOfRequestedColumnsPair = 0;

  // Populate outputColumns with columns that are requested from
  // the input dataset
  for (reqIt = this->Internals->Requests.begin();
    reqIt != this->Internals->Requests.end(); ++ reqIt)
    {
    // Each request contains only one pair of columns of interest
    // (if there are others, they are ignored).
    std::set<vtkStdString>::const_iterator colIt = reqIt->begin();
    const vtkStdString &colX = *colIt;
    if (!inData->GetColumnByName(colX.c_str()))
      {
      vtkWarningMacro("InData table does not have a column "
        << colX.c_str()
        << ". Ignoring this pair.");
      continue;
      }

    ++colIt;
    const vtkStdString &colY = *colIt;
    if (!inData->GetColumnByName(colY.c_str()))
      {
      vtkWarningMacro("InData table does not have a column "
        << colY.c_str()
        << ". Ignoring this pair.");
      continue;
      }

    // Verify column types
    vtkDataArray *inputColX =
      vtkDataArray::SafeDownCast(inData->GetColumnByName(colX.c_str()));
    vtkDataArray *inputColY =
      vtkDataArray::SafeDownCast(inData->GetColumnByName(colY.c_str()));
    if (!inputColX || !inputColY)
      {
      vtkErrorMacro(
        << "HDR cannot work with columns that are not of vtkDataArray type");
      return;
      }
    
    vtkSmartPointer<vtkDataArray> arrX =
      vtkDataArray::CreateDataArray(inputColX->GetDataType());
    arrX->DeepCopy(inputColX);
    arrX->SetName(inputColX->GetName());
    outputColumns->AddColumn(arrX);

    vtkSmartPointer<vtkDataArray> arrY =
      vtkDataArray::CreateDataArray(inputColY->GetDataType());
    arrY->DeepCopy(inputColY);
    arrY->SetName(inputColY->GetName());
    outputColumns->AddColumn(arrY);

    // Compute for the two columns and each observations the estimator of
    // density. Create a double Array that contains number of requested data
    // series components. Each tuple will contain the correspondent value
    // casted if necessary into a double.

    vtkNew<vtkDoubleArray> inObservations;
    inObservations->SetNumberOfComponents(2);
    inObservations->SetNumberOfTuples(outputColumns->GetNumberOfRows());

    inObservations->CopyComponent(0, inputColX, 0);
    inObservations->CopyComponent(1, inputColY, 0);

    // outObservations store the density vector
    vtkSmartPointer<vtkDataArray> outObservations =
      vtkDataArray::CreateDataArray(inObservations->GetDataType());
    outObservations->SetNumberOfComponents(1);
    outObservations->SetNumberOfTuples(inObservations->GetNumberOfTuples());

    this->ComputeHDR(inObservations.GetPointer(), outObservations);
    std::stringstream name;
    name << "HDR"; //f(" << inputColX->GetName() << "," << inputColY->GetName() << ")";
    outObservations->SetName(name.str().c_str());
    outputColumns->AddColumn(outObservations);

    // One requested pair of columns has been added.
    this->NumberOfRequestedColumnsPair += 1;
    } // End requests iteration.

  outMeta->SetNumberOfBlocks(1);
  outMeta->SetBlock(0, outputColumns.GetPointer());
  outMeta->GetMetaData(static_cast<unsigned int>(0))->
    Set(vtkCompositeDataSet::NAME(), "Estimator of density Data");
  }

// ----------------------------------------------------------------------
void vtkHighestDensityRegionsStatistics::Derive(vtkMultiBlockDataSet* outMeta)
  {
  vtkTable *hdrTab;

  if (!outMeta ||
    !(hdrTab = vtkTable::SafeDownCast(outMeta->GetBlock(0)))
    )
    {
    vtkErrorMacro(<< "Cannot derive: learn method never called");
    return;
    }

  // Compute sorted HDR for the last NumberOfRequestedColumnsPair
  // columns computed by Learn method. Sorted HDR is of same type as HDR
  for (vtkIdType i = 0; i < this->NumberOfRequestedColumnsPair; ++i)
    {
    vtkDataArray *hdr;
    if (!(hdr = vtkDataArray::SafeDownCast(
      hdrTab->GetColumn(hdrTab->GetNumberOfColumns()
      - this->NumberOfRequestedColumnsPair + i))))
      {
      vtkErrorMacro(<< "HDR column must be of vtkDataArray type");
      return;
      }

    double normFactor = 0.0;
    for (vtkIdType j = 0; j < hdr->GetNumberOfTuples(); ++j)
      {
      normFactor += hdr->GetTuple1(j);
      }

    if (normFactor != 0.0)
      {
      normFactor = 1.0 / normFactor;
      }

    // Creation of the hdr array.
    vtkSmartPointer<vtkDataArray> normalizedHDR =
      vtkDataArray::CreateDataArray(hdr->GetDataType());
    std::string name = "HDRn";//std::string(hdr->GetName()) + "N";
    normalizedHDR->SetName(name.c_str());
    normalizedHDR->SetNumberOfComponents(1);
    normalizedHDR->SetNumberOfTuples(hdr->GetNumberOfTuples());
    for (vtkIdType j = 0; j < hdr->GetNumberOfTuples(); ++j)
      {
      normalizedHDR->SetTuple1(j, normFactor * hdr->GetTuple1(j));
      }

    hdrTab->AddColumn(normalizedHDR);
    }
  }

// ----------------------------------------------------------------------
double vtkHighestDensityRegionsStatistics::ComputeHDR(vtkDataArray *inObservations,
                                                      vtkDataArray *outDensity)
  {
  vtkIdType nbObservations = inObservations->GetNumberOfTuples();

  if (nbObservations == 0)
    {
    vtkErrorMacro(<< "Empty observation array");
    return 0.0;
    }
  double sum  = 0.0;

  double denom = 1.0 / static_cast<double>(nbObservations);

  // Let's compute the HDR for each points of the observations
  for (vtkIdType i = 0; i < nbObservations; ++i)
    {
    double currentXi[2];
    double currentXj[2];
    double hdr = 0.0;

    // We are working in a bivariate model.
    inObservations->GetTuple(i, currentXi);
    // Sum all gaussian kernel
    for (vtkIdType j = 0; j < nbObservations; ++j)
      {
      // Avoid case where point is compared to itself
      if (i == j)
        {
        continue;
        }
      inObservations->GetTuple(j, currentXj);

      hdr += this->ComputeSmoothGaussianKernel(
        inObservations->GetNumberOfComponents(),
        currentXi[0] - currentXj[0],
        currentXi[1] - currentXj[1]);
      }
    double d = denom * hdr;
    outDensity->SetTuple1(i, d);
    sum += d;
    }

  return sum;
  }

// ----------------------------------------------------------------------
double vtkHighestDensityRegionsStatistics::ComputeSmoothGaussianKernel(
  int dimension, double khx, double khy)
  {
  double HDeterminant =
    vtkMath::Determinant2x2(this->SmoothHC1, this->SmoothHC2);
  // TODO: Check if HDeterminant is negative ! and what to do...
  if (HDeterminant != 0.0)
    {
    HDeterminant = 1.0 / sqrt(HDeterminant);
    }

  // We need to multiply the input vector by the smooth square root of
  // H matrix parameter: sqrt(H) * [khx, khy] -> random vector of the
  // standard gaussian input.

  // If a H coefficient is equal to 0.0. we don't compute its sqrt to avoid
  // domain error.
  double SHC10 = 0.0;
  double SHC11 = 0.0;
  double SHC20 = 0.0;
  double SHC21 = 0.0;

  if (this->SmoothHC1[0] != 0.0)
    {
    SHC10 = 1.0 / sqrt(this->SmoothHC1[0]);
    }
  if (this->SmoothHC1[1] != 0.0)
    {
    SHC11 = 1.0 / sqrt(this->SmoothHC1[1]);
    }
  if (this->SmoothHC2[0] != 0.0)
    {
    SHC20 = 1.0 / sqrt(this->SmoothHC2[0]);
    }
  if (this->SmoothHC2[1] != 0.0)
    {
    SHC21 = 1.0 / sqrt(this->SmoothHC2[1]);
    }

  // Call the standard gaussian kernel with the new random vector.
  return HDeterminant *
    this->ComputeStandardGaussianKernel(dimension,
    SHC10 * khx + SHC11 * khy,
    SHC20 * khx + SHC21 * khy);
  }

// ----------------------------------------------------------------------
double vtkHighestDensityRegionsStatistics::ComputeStandardGaussianKernel(
  int vtkNotUsed(dimension), double kx, double ky)
  {
  return exp(-(kx * kx + ky * ky) / 2.0) / (2.0 * vtkMath::Pi());
  }