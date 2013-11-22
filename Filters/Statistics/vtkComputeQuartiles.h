/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkComputeQuartiles.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __vtkComputeQuartiles_h
#define __vtkComputeQuartiles_h

#include "vtkFiltersStatisticsModule.h" // For export macro
#include "vtkTableAlgorithm.h"

class vtkDataSet;
class vtkDoubleArray;
class vtkFieldData;
class vtkTable;

// .NAME vtkComputeQuartiles - Extract quartiles from any dataset
// .SECTION Description
// vtkComputeQuartiles accepts any vtkDataSet as input and produces a
// vtkTable data as output.

class VTKFILTERSSTATISTICS_EXPORT vtkComputeQuartiles : public vtkTableAlgorithm
{
public:
  static vtkComputeQuartiles* New();
  vtkTypeMacro(vtkComputeQuartiles, vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
    
protected: 
  vtkComputeQuartiles();
  ~vtkComputeQuartiles();

  virtual int FillInputPortInformation (int port, vtkInformation *info);

  virtual int RequestData(vtkInformation *request, 
                          vtkInformationVector **inputVector, 
                          vtkInformationVector *outputVector); 

  void ComputeTable(vtkDataObject*, vtkTable*, vtkIdType);

  int FieldAssociation;
    
private:
  void operator=(const vtkComputeQuartiles&); // Not implemented
  vtkComputeQuartiles(const vtkComputeQuartiles&); // Not implemented
  
  int GetInputFieldAssociation();
  vtkFieldData* GetInputFieldData(vtkDataObject* input);
};

#endif
