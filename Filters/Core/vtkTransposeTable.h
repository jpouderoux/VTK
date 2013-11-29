/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTransposeTable.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkTransposeTable- Computes the transpose of an input table.
//
// .SECTION Thanks
// Developed by Joachim Pouderoux at Kitware SAS 2013

#ifndef __vtkTransposeTable_h
#define __vtkTransposeTable_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkTableAlgorithm.h"

class VTKFILTERSCORE_EXPORT vtkTransposeTable : public vtkTableAlgorithm
{
public:
  static vtkTransposeTable* New();
  vtkTypeMacro(vtkTransposeTable, vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // This flag indicates if a column must be inserted at index 0
  // with the names (ids) of the input columns. 
  // Default: true
  vtkGetMacro(AddIdColumn, bool);
  vtkSetMacro(AddIdColumn, bool);
  vtkBooleanMacro(AddIdColumn, bool);

  // Description:
  // This flag indicates if the output column must be named using the
  // names listed in the index 0 column.
  // Default: false
  vtkGetMacro(UseIdColumn, bool);
  vtkSetMacro(UseIdColumn, bool);
  vtkBooleanMacro(UseIdColumn, bool);
    
  // Description:
  // Get/Set the name of the id column added by option AddIdColumn.
  // Default: ColName
  vtkGetStringMacro(IdColumnName);
  vtkSetStringMacro(IdColumnName);

protected:
  vtkTransposeTable();
  ~vtkTransposeTable();

  int RequestData(vtkInformation*,
    vtkInformationVector**,
    vtkInformationVector*);

  bool AddIdColumn;
  bool UseIdColumn;
  char* IdColumnName;

private:
  vtkTransposeTable(const vtkTransposeTable&); // Not implemented
  void operator=(const vtkTransposeTable&);   // Not implemented
};

#endif
