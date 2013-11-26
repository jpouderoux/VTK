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
#include "vtkChartLegend.h"
#include "vtkContextScene.h"
#include "vtkContextView.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkLookupTable.h"
#include "vtkNew.h"
#include "vtkPlotBox.h"
#include "vtkRect.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTable.h"

//----------------------------------------------------------------------------
int TestBoxPlot(int , char * [])
{
  // Set up a 2D scene, add an XY chart to it
  vtkNew<vtkContextView> view;
  view->GetRenderWindow()->SetSize(400, 400);
  view->GetRenderWindow()->SetMultiSamples(0);

  vtkNew<vtkChartXY> chart;
  view->GetScene()->AddItem(chart.GetPointer());
  chart->GetAxis(vtkAxis::BOTTOM)->SetLabelsVisible(false);
  chart->GetAxis(vtkAxis::BOTTOM)->SetTicksVisible(false);  
  chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("");  
  chart->GetAxis(vtkAxis::LEFT)->SetTitle("");
  chart->GetLegend()->SetVerticalAlignment(vtkChartLegend::BOTTOM);
  chart->SetShowLegend(true);

  // Creates a vtkPlotBox input table
  // The vtkPlotBox object will display 4 (arbitrary) box plot
  int numParam = 4;
  vtkNew<vtkTable> inputBoxPlotTable;
  
  for (int i = 0; i < numParam; i++)
    {
    char num[3];
    sprintf(num, "%d", i);
    char name[10];
    strcpy(name,"Param ");
    strcat(name,num);
  
    vtkNew<vtkIntArray> arrIndex;
    arrIndex->SetName(name);
    inputBoxPlotTable->AddColumn(arrIndex.GetPointer());
    }

  inputBoxPlotTable->SetNumberOfRows(5);
    
  for (int i = 0; i < numParam; i++)
    {    
    inputBoxPlotTable->SetValue(0, i, i/2);     //Q0
    inputBoxPlotTable->SetValue(1, i, 2*i + 2 - i); //Q1
    inputBoxPlotTable->SetValue(2, i, 2*i + 4); //Q2
    inputBoxPlotTable->SetValue(3, i, 2*i + 7); //Q3
    inputBoxPlotTable->SetValue(4, i, 2*i + 8); //Q4
    }

  vtkNew<vtkLookupTable> lookup;
  lookup->SetNumberOfColors(5);
  lookup->SetRange(0, 4);
  lookup->Build();   
  
  for (int i = 0; i < numParam; i++)
    {
    // Set the vtkPlotBox (here 4 box have to been displayed)
    vtkNew<vtkPlotBox> boxPlot;
    chart->AddPlot(boxPlot.GetPointer());
    boxPlot->SetUseIndexForXSeries(true);
    boxPlot->SetInputData(inputBoxPlotTable.GetPointer(), 0, i);
    double rgb[3];
    lookup->GetColor(i, rgb);
    boxPlot->SetColor(rgb[0], rgb[1], rgb[2]);
    }
  
  // Render the scene
  view->GetRenderWindow()->SetMultiSamples(0);
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();

  return EXIT_SUCCESS;
}
