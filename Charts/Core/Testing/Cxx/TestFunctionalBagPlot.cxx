/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestLinePlot.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkChartXY.h"
#include "vtkContextScene.h"
#include "vtkContextView.h"
#include "vtkDoubleArray.h"
#include "vtkPlotFunctionalBag.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStringArray.h"
#include "vtkTable.h"

#include <sstream>

//----------------------------------------------------------------------------
int TestFunctionalBagPlot(int, char * [])
{
  // Set up a 2D scene, add an XY chart to it
  vtkNew<vtkContextView> view;
  view->GetRenderWindow()->SetSize(400, 400);
  view->GetRenderWindow()->SetMultiSamples(0);
  vtkNew<vtkChartXY> chart;
  view->GetScene()->AddItem(chart.GetPointer());
  chart->SetShowLegend(true);

  // Creates an input table  
  const int numCols = 7;
  const int numVals = 100; 

  vtkNew<vtkTable> inputTable;
  for (int i = 0; i < numCols; i++)
    {
    vtkNew<vtkDoubleArray> arr;
    std::stringstream ss;
    ss << "Y" << i;
    arr->SetName(ss.str().c_str());
    arr->SetNumberOfValues(numVals);
    inputTable->AddColumn(arr.GetPointer());
    for (int j = 0; j < numVals; j++)
      {
      arr->SetValue(j, (i+1) * abs(sin(j*(2*vtkMath::Pi())/(float)numVals)) * j + i*20);//rand()/(double)RAND_MAX);  
      }
    }

  vtkNew<vtkTable> inputDensityTable;
  vtkNew<vtkDoubleArray> arr;
  arr->SetName("Density");
  arr->SetNumberOfValues(numCols);
  arr->SetValue(0, 0.08);
  arr->SetValue(1, 0.08);
  arr->SetValue(2, 0.12);
  arr->SetValue(3, 0.25);
  arr->SetValue(4, 0.25);
  arr->SetValue(5, 0.12);
  arr->SetValue(6, 0.08);

  inputDensityTable->AddColumn(arr.GetPointer());
  
  vtkNew<vtkStringArray> varArr;
  varArr->SetName("Variable");
  varArr->SetNumberOfValues(numCols);
  for (int j = 0; j < numCols; j++)
    {
    std::stringstream ss;
    ss << "Y" << j;
    varArr->SetValue(j, ss.str().c_str());
    }
  inputDensityTable->AddColumn(varArr.GetPointer());

  // Create the functional bag plot
  vtkNew<vtkPlotFunctionalBag> plot;
  chart->AddPlot(plot.GetPointer());
  plot->SetInputData(inputTable.GetPointer());
  plot->SetInputDensityData(inputDensityTable.GetPointer(),"Density", "Variable");  
  plot->SetColor(255, 0, 0, 255);
  
  // Render the scene  
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();

  return EXIT_SUCCESS;
}
