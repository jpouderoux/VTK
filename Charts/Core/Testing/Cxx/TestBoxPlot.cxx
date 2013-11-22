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

#include "vtkAxis.h"
#include "vtkChartXY.h"
#include "vtkContextScene.h"
#include "vtkContextView.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkNew.h"
#include "vtkOrderStatistics.h"
#include "vtkPlotBox.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTable.h"
#include "vtkTimerLog.h"

//----------------------------------------------------------------------------
int TestBoxPlot(int , char * [])
{
  // Set up a 2D scene, add an XY chart to it
  vtkNew<vtkContextView> view;
  view->GetRenderWindow()->SetSize(800, 600);
  vtkNew<vtkChartXY> chart;
  view->GetScene()->AddItem(chart.GetPointer());

  // Creates a vtkPlotBox input table
  // The vtkPlotBox object will display 4 (arbitrary) box plot
  int numData = 4;
  vtkNew<vtkTable> inputBoxPlotTable;

  vtkNew<vtkIntArray> arrIndex;
  arrIndex->SetName("Index");
  inputBoxPlotTable->AddColumn(arrIndex.GetPointer());

  vtkNew<vtkDoubleArray> arrQ0;
  arrQ0->SetName("Q0");
  inputBoxPlotTable->AddColumn(arrQ0.GetPointer());

  vtkNew<vtkDoubleArray> arrQ1;
  arrQ1->SetName("Q1");
  inputBoxPlotTable->AddColumn(arrQ1.GetPointer());

  vtkNew<vtkDoubleArray> arrQ2;
  arrQ2->SetName("Q2");
  inputBoxPlotTable->AddColumn(arrQ2.GetPointer());

  vtkNew<vtkDoubleArray> arrQ3;
  arrQ3->SetName("Q3");
  inputBoxPlotTable->AddColumn(arrQ3.GetPointer());

  vtkNew<vtkDoubleArray> arrQ4;
  arrQ4->SetName("Q4");
  inputBoxPlotTable->AddColumn(arrQ4.GetPointer());

  // Test charting with numData box plot.
  inputBoxPlotTable->SetNumberOfRows(numData);

  // Fill inputBoxPlotTable
  // inputBoxPlotTable need at least six column to works
  // One for the index representing x position and the others
  // for quartile values.
  for (int i = 0; i < numData; ++i)
    {
    inputBoxPlotTable->SetValue(i, 0, i); //C0
    inputBoxPlotTable->SetValue(i, 1, i); //C1
    inputBoxPlotTable->SetValue(i, 2, i + 2); //C2
    inputBoxPlotTable->SetValue(i, 3, i + 4); //C3
    inputBoxPlotTable->SetValue(i, 4, i + 7); //C4
    inputBoxPlotTable->SetValue(i, 5, i + 8); //C5
    }

  chart->SetShowLegend(true);

  // Set the vtkPlotBox (here 4 box have to been displayed)
  vtkNew<vtkPlotBox> boxPlot;
  chart->AddPlot(boxPlot.GetPointer());
  boxPlot->SetInputData(inputBoxPlotTable.GetPointer(), arrIndex->GetName(),
                        arrQ0->GetName(),
                        arrQ1->GetName(),
                        arrQ2->GetName(),
                        arrQ3->GetName(),
                        arrQ4->GetName());
  boxPlot->SetColor(255, 0, 0, 255);

  // Render the scene
  view->GetRenderWindow()->SetMultiSamples(0);
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();

  return EXIT_SUCCESS;
}