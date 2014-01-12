/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRGridSlice.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkRGridSlice.h"

#include "vtkBitArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataSetAttributes.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkRGrid.h"

vtkStandardNewMacro(vtkRGridSlice);

//-----------------------------------------------------------------------------
vtkRGridSlice::vtkRGridSlice()
{
}

//-----------------------------------------------------------------------------
vtkRGridSlice::~vtkRGridSlice()
{
}

//----------------------------------------------------------------------------
void vtkRGridSlice::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//-----------------------------------------------------------------------------
int vtkRGridSlice::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRGrid");
  return 1;
}

//----------------------------------------------------------------------------
int vtkRGridSlice::RequestData(vtkInformation*,
                               vtkInformationVector** inputVector,
                               vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Retrieve input and output
  vtkRGrid* input = vtkRGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkRGrid* output = vtkRGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int dims[3];
  input->GetDimensions(dims);

  // Initialize output cell data
  vtkDataSetAttributes* inData = static_cast<vtkDataSetAttributes*>(input->GetCellData());
  vtkDataSetAttributes* outData = static_cast<vtkDataSetAttributes*>(output->GetCellData());
  outData->CopyAllocate(inData);


  int mini = dims[0] / 2;
  int maxi = mini+1;

  output->SetPoints(input->GetPoints());
  output->SetDimensions(maxi - mini, dims[1], dims[2]);
  vtkNew<vtkCellArray> cells;
  cells->Allocate((maxi - mini) * dims[1] * dims[2] * 9);
  for (int k = 0; k < dims[2]; k++)
    {
    for (int j = 0; j < dims[1]; j++)
      {
      for (int i = mini; i < maxi; i++)
        {
        vtkIdType cellId = i + j * dims[0] + k * dims[0] * dims[1];
        vtkNew<vtkIdList> ptIds;
        input->GetCellPoints(cellId, ptIds.GetPointer());
        vtkIdType nCellId = cells->InsertNextCell(ptIds.GetPointer());
        outData->CopyData(inData, cellId, nCellId);
        }
      }
    }
  output->SetCells(cells.GetPointer());

  this->UpdateProgress(1.);

  return 1;
}
