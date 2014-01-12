/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRGridSlice.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkRGridSlice -
//
// .SECTION Description
//
// .SECTION See Also
// vtkRGrid
//
// .SECTION Thanks

#ifndef __vtkRGridSlice_h
#define __vtkRGridSlice_h

#include "vtkFiltersGeometryModule.h" // For export macro
#include "vtkRGridAlgorithm.h"

class vtkCellArray;
class vtkDataSetAttributes;
class vtkHyperTreeGrid;
class vtkPoints;

class VTKFILTERSGEOMETRY_EXPORT vtkRGridSlice : public vtkRGridAlgorithm
{
public:
  static vtkRGridSlice* New();
  vtkTypeMacro(vtkRGridSlice, vtkRGridAlgorithm);
  void PrintSelf(ostream&, vtkIndent);

  vtkGetVector6Macro(SliceExtents, int);
  vtkSetVector6Macro(SliceExtents, int);

protected:
  vtkRGridSlice();
  ~vtkRGridSlice();

  int SliceExtents[6];

  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  virtual int FillInputPortInformation(int, vtkInformation*);

private:
  vtkRGridSlice(const vtkRGridSlice&);  // Not implemented.
  void operator=(const vtkRGridSlice&);  // Not implemented.
};

#endif
