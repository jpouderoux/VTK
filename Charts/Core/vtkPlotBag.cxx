/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlotBag.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkBrush.h"
#include "vtkContext2D.h"
#include "vtkContextMapper2D.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkObjectFactory.h"
#include "vtkPen.h"
#include "vtkPlotBag.h"
#include "vtkPoints.h"
#include "vtkPoints2D.h"
#include "vtkPointsProjectedHull.h"
#include "vtkTable.h"
#include "vtkTimeStamp.h"

#include <algorithm>

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPlotBag);

//-----------------------------------------------------------------------------
vtkPlotBag::vtkPlotBag()
{
  this->MedianPoints = vtkPoints2D::New();
  this->Q1Points = vtkPoints2D::New();
}

//-----------------------------------------------------------------------------
vtkPlotBag::~vtkPlotBag()
{
  if (this->MedianPoints)
    {
    this->MedianPoints->Delete();
    this->MedianPoints = 0;
    }
  if (this->Q1Points)
    {
    this->Q1Points->Delete();
    this->Q1Points = 0;
    }
}

//-----------------------------------------------------------------------------
void vtkPlotBag::Update()
{
  this->Superclass::Update();

  if (!this->Visible)
    {
    return;
    }

  // Check if we have an input
  vtkTable *table = this->Data->GetInput();
  vtkDataArray *bag = vtkDataArray::SafeDownCast(
    this->Data->GetInputAbstractArrayToProcess(2, this->GetInput()));
  if (!table || !bag)
    {
    vtkDebugMacro(<< "Update event called with no input table or bag column set.");
    return;
    }
  if (this->Data->GetMTime() > this->BuildTime ||
    table->GetMTime() > this->BuildTime ||
    this->MTime > this->BuildTime)
    {
    vtkDebugMacro(<< "Updating cached values.");
    this->UpdateTableCache(bag);
    }  
}

//-----------------------------------------------------------------------------
class ArraySorter
{
public:
  ArraySorter(vtkDataArray* arr) : Array(arr) {}
  bool operator()(const vtkIdType& a, const vtkIdType& b)
  {
    return this->Array->GetTuple1(a) > this->Array->GetTuple1(b);
  }
  vtkDataArray* Array;
};

//-----------------------------------------------------------------------------
void vtkPlotBag::UpdateTableCache(vtkDataArray* density)
{
  vtkIdType nbPoints = density->GetNumberOfTuples();

  // Sort the density array
  std::vector<vtkIdType> ids;
  ids.resize(nbPoints);
  double sum = 0.0;
  for (vtkIdType i = 0; i < nbPoints; i++)
    {
    sum += density->GetTuple1(i);
    ids[i] = i;
    }
  
  vtkDoubleArray* nDensity = 0;
  // Normalize the density array if it is not
  if (fabs(sum - 1.0) > 1.0e-12)
    {
    sum = 1.0 / sum;
    nDensity = vtkDoubleArray::New();
    nDensity->SetNumberOfComponents(1);
    nDensity->SetNumberOfTuples(nbPoints);
    for (vtkIdType i = 0; i < nbPoints; i++)
      {
      nDensity->SetTuple1(i, density->GetTuple1(ids[i]) * sum);
      }
    density = nDensity;
    }
  
  // Sort array by density
  ArraySorter arraySorter(density);
  std::sort(ids.begin(), ids.end(), arraySorter);
  
  vtkNew<vtkPointsProjectedHull> q1Points;
  q1Points->Allocate(nbPoints);
  vtkNew<vtkPointsProjectedHull> medianPoints;
  medianPoints->Allocate(nbPoints);

  sum = 0.0;
  for (vtkIdType i = 0; i < nbPoints; i++)
    {
    sum += density->GetTuple1(ids[i]);
    if (sum > 0.75)
      {
      break;
      }
    //vtkErrorMacro( << "sum " << sum << " pt " << ids[i] );
    double x[3];
    this->Points->GetPoint(ids[i], x);
    if (sum <= 0.5)
      {
      medianPoints->InsertNextPoint(x);      
      //vtkErrorMacro( << "added in median" );
      }
    q1Points->InsertNextPoint(x);
    //vtkErrorMacro( << "added in q1" );
    }

  this->MedianPoints->Reset();
  this->Q1Points->Reset();

  // Compute the convex hull for the median points
  if (medianPoints->GetNumberOfPoints() > 0)
    {
    int size = medianPoints->GetSizeCCWHullZ();
    this->MedianPoints->SetDataTypeToFloat();
    this->MedianPoints->SetNumberOfPoints(size+1);
    medianPoints->GetCCWHullZ(
      (float*)this->MedianPoints->GetData()->GetVoidPointer(0), size);
    double x[3];
    this->MedianPoints->GetPoint(0, x);
    this->MedianPoints->SetPoint(size, x);
    }
 
  // Compute the convex hull for the first quartile points
  if (q1Points->GetNumberOfPoints() > 0)
    {
    int size = q1Points->GetSizeCCWHullZ();    
    this->Q1Points->SetDataTypeToFloat();
    this->Q1Points->SetNumberOfPoints(size+1);
    q1Points->GetCCWHullZ(
      (float*)this->Q1Points->GetData()->GetVoidPointer(0), size);
    double x[3];
    this->Q1Points->GetPoint(0, x);
    this->Q1Points->SetPoint(size, x);
    }

  if (nDensity)
    {
    nDensity->Delete();
    }
  this->BuildTime.Modified();
}

//-----------------------------------------------------------------------------
bool vtkPlotBag::Paint(vtkContext2D *painter)
{
  vtkDebugMacro(<< "Paint event called in vtkPlotBag.");

  vtkTable *table = this->Data->GetInput();
  vtkDataArray *density = vtkDataArray::SafeDownCast(
    this->Data->GetInputAbstractArrayToProcess(2, this->GetInput()));
  if (!this->Visible || !this->Points || !table)
    {
    return false;
    }
    
  if (density)
    {
    // Let's draw the bags
    unsigned char pcolor[4];
    this->Pen->GetColor(pcolor);
    this->Pen->SetColor(0, 0, 0);
    unsigned char bcolor[4];
    this->Brush->GetColor(bcolor);
    this->Brush->SetColor(bcolor[0] / 2, bcolor[1] / 2, bcolor[2] / 2);
    painter->ApplyPen(this->Pen);
    painter->ApplyBrush(this->Brush);
    if (this->Q1Points->GetNumberOfPoints() > 0)
      {
      painter->DrawPolygon(this->Q1Points);
      }     
    this->Brush->SetColor(bcolor);
    
    painter->ApplyPen(this->Pen);
    painter->ApplyBrush(this->Brush);
    this->Brush->SetOpacity(128);
    if (this->MedianPoints->GetNumberOfPoints() > 0)
      {
      painter->DrawPolygon(this->MedianPoints);
      }

    this->Pen->SetColor(pcolor);
    }
  
  // Let PlotPoints draw the points
  return this->Superclass::Paint(painter);
}

//-----------------------------------------------------------------------------
bool vtkPlotBag::PaintLegend(vtkContext2D *painter, const vtkRectf& rect, int)
{
  painter->ApplyPen(this->Pen);
  painter->DrawLine(
    rect[0], rect[1] + 0.5 * rect[3],
    rect[0] + rect[2], rect[1] + 0.5 * rect[3]);

  return this->Superclass::PaintLegend(painter, rect, 0);
}

//-----------------------------------------------------------------------------
void vtkPlotBag::SetInputData(vtkTable *table)
{
  this->Data->SetInputData(table);
  this->Modified();
}

//-----------------------------------------------------------------------------
void vtkPlotBag::SetInputData(vtkTable *table, const vtkStdString &yColumn,
                              const vtkStdString &densityColumn)
{

  vtkDebugMacro(<< "Setting input, Y column = \"" << yColumn.c_str() << "\", "
                << "Density column = \"" << densityColumn.c_str() << "\"");

  if (table->GetColumnByName(densityColumn.c_str())->GetNumberOfTuples() 
    != table->GetColumnByName(yColumn.c_str())->GetNumberOfTuples())
    {
    vtkErrorMacro(<< "Input table not correctly initialized!");
    return;
    }

  this->SetInputData(table, yColumn, yColumn, densityColumn);
  this->UseIndexForXSeries = true;
}

//-----------------------------------------------------------------------------
void vtkPlotBag::SetInputData(vtkTable *table, const vtkStdString &xColumn,
                              const vtkStdString &yColumn,
                              const vtkStdString &densityColumn)
{
  vtkDebugMacro(<< "Setting input, X column = \"" << xColumn.c_str()
                << "\", " << "Y column = \"" 
                << yColumn.c_str() << "\""
                << "\", " << "Density column = \"" 
                << densityColumn.c_str() << "\"");

  this->Data->SetInputData(table);  
  this->Data->SetInputArrayToProcess(0, 0, 0,
    vtkDataObject::FIELD_ASSOCIATION_ROWS, xColumn.c_str());
  this->Data->SetInputArrayToProcess(1, 0, 0,
    vtkDataObject::FIELD_ASSOCIATION_ROWS, yColumn.c_str());
  this->Data->SetInputArrayToProcess(2, 0, 0,
    vtkDataObject::FIELD_ASSOCIATION_ROWS, densityColumn.c_str());
}

//-----------------------------------------------------------------------------
void vtkPlotBag::SetInputData(vtkTable *table, vtkIdType xColumn,
                              vtkIdType yColumn,
                              vtkIdType densityColumn)
{
  this->SetInputData(table,
    table->GetColumnName(xColumn), 
    table->GetColumnName(yColumn),
    table->GetColumnName(densityColumn));
}

//-----------------------------------------------------------------------------
void vtkPlotBag::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
