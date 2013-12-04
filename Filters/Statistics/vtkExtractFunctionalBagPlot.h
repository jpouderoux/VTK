/*=========================================================================

  Program:   ParaView
  Module:    vtkExtractFunctionalBagPlot.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkExtractFunctionalBagPlot
//
// .SECTION Description

#ifndef __vtkExtractFunctionalBagPlot_h
#define __vtkExtractFunctionalBagPlot_h

#include "vtkFiltersStatisticsModule.h" // For export macro
#include "vtkTableAlgorithm.h"


class VTKFILTERSSTATISTICS_EXPORT vtkExtractFunctionalBagPlot : public vtkTableAlgorithm
{
public:
  static vtkExtractFunctionalBagPlot* New();
  vtkTypeMacro(vtkExtractFunctionalBagPlot, vtkTableAlgorithm);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkExtractFunctionalBagPlot();
  virtual ~vtkExtractFunctionalBagPlot();

  int RequestData(vtkInformation*,
    vtkInformationVector**,
    vtkInformationVector*);

private:
  vtkExtractFunctionalBagPlot( const vtkExtractFunctionalBagPlot& ); // Not implemented.
  void operator = ( const vtkExtractFunctionalBagPlot& ); // Not implemented.
};

#endif // __vtkExtractFunctionalBagPlot_h
