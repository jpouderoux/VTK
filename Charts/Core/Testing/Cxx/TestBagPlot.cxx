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
#include "vtkIntArray.h"
#include "vtkNew.h"
#include "vtkPlotBag.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTable.h"

//----------------------------------------------------------------------------
int TestBagPlot(int, char * [])
{
  // Set up a 2D scene, add an XY chart to it
  vtkNew<vtkContextView> view;
  view->GetRenderWindow()->SetSize(400, 400);
  view->GetRenderWindow()->SetMultiSamples(0);
  vtkNew<vtkChartXY> chart;
  view->GetScene()->AddItem(chart.GetPointer());

  // Creates a vtkPlotBag input table
  // We construct a 2D grid 20*20.
  // the bag will represent a square inside of side 5
  int numDataI = 20;
  int numDataJ = 20;

  vtkNew<vtkTable> inputBagPlotTable;

  vtkNew<vtkIntArray> arrX;
  arrX->SetName("X");
  inputBagPlotTable->AddColumn(arrX.GetPointer());

  vtkNew<vtkDoubleArray> arrY;
  arrY->SetName("Y");
  inputBagPlotTable->AddColumn(arrY.GetPointer());

  vtkNew<vtkDoubleArray> arrDensity;
  arrDensity->SetName("Density");
  inputBagPlotTable->AddColumn(arrDensity.GetPointer());

  // Test charting with numData Bag plot.
  inputBagPlotTable->SetNumberOfRows(numDataI * numDataJ);

  // Fill inputBagPlotTable
  // inputBagPlotTable need at least six columns to works
  // One for the index representing x position and the others
  // for quartile values.

  for (int j = 0; j < numDataJ; ++j)
    {
    for (int i = 0; i < numDataI; ++i)
      {
      inputBagPlotTable->SetValue(j * numDataI + i, 0, i); //X
      inputBagPlotTable->SetValue(j * numDataI + i, 1, j); //Y
      inputBagPlotTable->SetValue(j * numDataI + i, 2, rand()/(double)RAND_MAX); // Density
      }
    }
  
  chart->SetShowLegend(true);

  vtkNew<vtkPlotBag> bagPlot;
  chart->AddPlot(bagPlot.GetPointer());
  bagPlot->SetInputData(inputBagPlotTable.GetPointer(), arrX->GetName(),
    arrY->GetName(), arrDensity->GetName());
  bagPlot->SetColor(255, 0, 0, 255);

  // Render the scene  
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();

  return EXIT_SUCCESS;
}
