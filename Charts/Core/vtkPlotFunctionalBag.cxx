/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlotFunctionalBag.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPlotFunctionalBag.h"

#include "vtkAbstractArray.h"
#include "vtkAxis.h"
#include "vtkBrush.h"
#include "vtkContext2D.h"
#include "vtkContextMapper2D.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkPen.h"
#include "vtkPlot.h"
#include "vtkPlotLine.h"
#include "vtkPoints2D.h"
#include "vtkRect.h"
#include "vtkScalarsToColors.h"
#include "vtkStringArray.h"
#include "vtkTable.h"

#include <algorithm>
#include <map>
#include <vector>

//-----------------------------------------------------------------------------
class vtkPlotFuntionalBagInternal
{
public:
  vtkPlotFuntionalBagInternal(vtkPlotFunctionalBag* parent) : Parent(parent) {}
  ~vtkPlotFuntionalBagInternal();
  
  std::vector<vtkPlotLine*> Lines;

protected:
  vtkPlotFunctionalBag* Parent;
};


//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPlotFunctionalBag);

//-----------------------------------------------------------------------------
vtkPlotFunctionalBag::vtkPlotFunctionalBag()
{
  this->LogX = false;
  this->LogY = false;

  this->Internal = new vtkPlotFuntionalBagInternal(this);
  this->LookupTable = 0;

  this->MedianPoints = vtkPoints2D::New();
  this->Q3Points = vtkPoints2D::New();
}

//-----------------------------------------------------------------------------
vtkPlotFunctionalBag::~vtkPlotFunctionalBag()
{
  if (this->LookupTable)
    {
    this->LookupTable->UnRegister(this);
    }

  if (this->MedianPoints)
    {
    this->MedianPoints->Delete();
    }
  if (this->Q3Points)
    {
    this->Q3Points->Delete();
    }
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::SetInputData(vtkTable *table)
{
  this->Data->SetInputData(table);  
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::SetInputDensityData(vtkTable *table, 
                                               const vtkStdString &densityColumn,
                                               const vtkStdString &variableNameColumn)
{
  this->Data->SetInputDataObject(1, table);
  this->Data->SetInputArrayToProcess(3, 1, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     densityColumn.c_str());
  if (variableNameColumn != "")
    {
    this->Data->SetInputArrayToProcess(4, 1, 0,
                                       vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                       variableNameColumn.c_str());
    }
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::SetInputDensityData(vtkTable *table, 
                                               vtkIdType densityColumn,
                                               vtkIdType variableNameColumn)
{
  this->SetInputDensityData(table,
                            table->GetColumnName(densityColumn),
                            variableNameColumn >= 0 ? 
                            table->GetColumnName(variableNameColumn) : "");
}

//-----------------------------------------------------------------------------
class DensityVal
{
public:
  DensityVal(double d, vtkAbstractArray* arr) : Density(d), Array(arr) {}
  bool operator<(const DensityVal& b)
  {
    return this->Density > b.Density;
  }
  double Density;
  vtkAbstractArray* Array;
};

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::Update()
{
  if (!this->Visible)
    {
    return;
    }
  // Check if we have an input
  vtkTable *table = this->Data->GetInput();

  if (!table)
    {
    vtkDebugMacro(<< "Update event called with no input table set.");
    return;
    }
  else if(this->Data->GetMTime() > this->BuildTime ||
          table->GetMTime() > this->BuildTime ||
          (this->LookupTable && this->LookupTable->GetMTime() > this->BuildTime) ||
          this->MTime > this->BuildTime)
    {
    vtkDebugMacro(<< "Updating cached values.");
    this->UpdateTableCache(table);

    vtkTable *densityTable = 
      vtkTable::SafeDownCast(this->Data->GetInputDataObject(1, 0));
    if (!densityTable)
      {
      vtkDebugMacro(<< "Update event called with no input density table set.");
      return;
      }
    vtkDoubleArray *density = vtkDoubleArray::SafeDownCast(
      this->Data->GetInputAbstractArrayToProcess(3, densityTable));
    if (!density)
      {
      vtkDebugMacro(<< "Update event called with non double density array.");
      return;
      }

    vtkStringArray *varName = vtkStringArray::SafeDownCast(
      this->Data->GetInputAbstractArrayToProcess(4, densityTable));
    if (!varName)
      {
      vtkDebugMacro(<< "Update event called with non double density array.");
      return;
      }

    // Fetch and sort arrays according their density
    std::vector<DensityVal> varNames;
    for (int i = 0; i < varName->GetNumberOfValues(); i++)
      {
      varNames.push_back(
        DensityVal(density->GetValue(i), table->GetColumnByName(varName->GetValue(i))));
      }

    std::sort(varNames.begin(), varNames.end());
    
    std::vector<vtkAbstractArray*> medianLines;
    std::vector<vtkAbstractArray*> q3Lines;
    double sum = 0.0;
    for (size_t i = 0; i < varNames.size(); i++)
      {
      sum += varNames[i].Density;
      if (sum <= 0.75)
        {
        if (sum <= 0.5)
          {
          medianLines.push_back(varNames[i].Array);
          }
        q3Lines.push_back(varNames[i].Array);
        }
      }

    // Generate the quad strip arrays
    vtkIdType nbRows = table->GetNumberOfRows();
    this->MedianPoints->Reset();
    this->MedianPoints->SetNumberOfPoints(2 * nbRows + 0);
    this->Q3Points->Reset();
    this->Q3Points->SetNumberOfPoints(2 * nbRows + 0);

    size_t medianCount = medianLines.size();
    size_t q3Count = q3Lines.size();
    for (vtkIdType i = 0; i < nbRows; i++)
      {
      double vMin = VTK_DOUBLE_MAX;
      double vMax = VTK_DOUBLE_MIN;
      for (size_t j = 0; j < medianCount; j++)
        {
        double v = medianLines[j]->GetVariantValue(i).ToDouble();
        if (v < vMin) { vMin = v; }
        if (v > vMax) { vMax = v; }
        }
      this->MedianPoints->SetPoint(2 * i, i, vMin);
      this->MedianPoints->SetPoint(2 * i + 1, i, vMax);

      vMin = VTK_DOUBLE_MAX;
      vMax = VTK_DOUBLE_MIN;
      for (size_t j = 0; j < q3Count; j++)
        {
        double v = q3Lines[j]->GetVariantValue(i).ToDouble();
        if (v < vMin) { vMin = v; }
        if (v > vMax) { vMax = v; }
        }
      this->Q3Points->SetPoint(2*i, i, vMin);
      this->Q3Points->SetPoint(2*i + 1, i, vMax);
      }
    }
}

//-----------------------------------------------------------------------------
bool vtkPlotFunctionalBag::UpdateTableCache(vtkTable *table)
{
  for (unsigned int i = 0; i < this->Internal->Lines.size(); i++)
    {
    this->RemoveItem(this->Internal->Lines[i]);
    }
  this->Internal->Lines.clear();

  if (!this->LookupTable)
    {
    this->CreateDefaultLookupTable();
    this->LookupTable->SetRange(0, table->GetNumberOfColumns());
    this->LookupTable->Build();
    }

  for (int i = 0; i < table->GetNumberOfColumns(); i++)
    {
    vtkNew<vtkPlotLine> line;
    line->SetInputData(table, 0, i);
    line->SetUseIndexForXSeries(true);
    line->SetMarkerStyle(vtkPlotPoints::NONE);
    double rgb[3];
    this->LookupTable->GetColor(i, rgb);
    line->SetColor(rgb[0], rgb[1], rgb[2]);
    line->SetWidth(this->GetWidth());    
    this->AddItem(line.GetPointer());
    this->Internal->Lines.push_back(line.GetPointer());
    }

  for (unsigned int i = 0; i < this->GetNumberOfItems(); i++)
    {
    vtkPlot::SafeDownCast(this->GetItem(i))->Update();
    }
  this->BuildTime.Modified();

  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotFunctionalBag::Paint(vtkContext2D *painter)
{
  // This is where everything should be drawn, or dispatched to other methods.
  vtkDebugMacro(<< "Paint event called in vtkPlotFunctionalBag.");

  if (!this->Visible)
    {
    return false;
    }

  // Let's draw the bags
  unsigned char pcolor[4];
  double pwidth = this->Pen->GetWidth();
  this->Pen->SetWidth(9);
  this->Pen->GetColor(pcolor);
  this->Pen->SetColor(0, 0, 0);
  unsigned char bcolor[4];
  this->Brush->SetColor(pcolor);
  this->Brush->GetColor(bcolor);
  this->Brush->SetColor(bcolor[0] / 2, bcolor[1] / 2, bcolor[2] / 2);
  painter->ApplyPen(this->Pen);
  painter->ApplyBrush(this->Brush);
  if (this->Q3Points->GetNumberOfPoints() > 0)
    {
    painter->DrawQuadStrip(this->Q3Points);
    }

  this->Brush->SetColor(bcolor);
  this->Pen->SetColor(pcolor);

  painter->ApplyPen(this->Pen);
  painter->ApplyBrush(this->Brush);
  this->Brush->SetOpacity(128);
  if (this->MedianPoints->GetNumberOfPoints() > 0)
    {
    painter->DrawQuadStrip(this->MedianPoints);    
    }

  this->Pen->SetWidth(pwidth);

  this->PaintChildren(painter);

  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotFunctionalBag::PaintLegend(vtkContext2D *painter, 
                                       const vtkRectf& rect, int index)
{
  this->Internal->Lines[index]->PaintLegend(painter, rect, index);
  return true;
}

//-----------------------------------------------------------------------------
vtkStringArray* vtkPlotFunctionalBag::GetLabels()
{
  // If the label string is empty, return the y column name
  if (this->Labels)
    {
    return this->Labels;
    }
  else if (this->AutoLabels)
    {
    return this->AutoLabels;
    }
  else if (this->Data->GetInput())
    {
    this->AutoLabels = vtkSmartPointer<vtkStringArray>::New();
    for (int i = 0; i < this->Internal->Lines.size(); i++)
      {
      this->AutoLabels->InsertNextValue(this->Internal->Lines[i]->GetLabel());
      }
    return this->AutoLabels;
    }
  else
    {
    return NULL;
    }
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::GetBounds(double bounds[4])
{
  bounds[0] = VTK_DOUBLE_MAX;
  bounds[1] = VTK_DOUBLE_MIN;
  bounds[2] = VTK_DOUBLE_MAX;
  bounds[3] = VTK_DOUBLE_MIN;
  for (unsigned int i = 0; i < this->Internal->Lines.size(); i++)
    {
    double b[4];
    this->Internal->Lines[i]->GetBounds(b);
    if (b[0] < bounds[0]) { bounds[0] = b[0]; }
    if (b[1] > bounds[1]) { bounds[1] = b[1]; }
    if (b[2] < bounds[2]) { bounds[2] = b[2]; }
    if (b[3] > bounds[3]) { bounds[3] = b[3]; }
    }
  vtkDebugMacro(<< "Bounds: " << bounds[0] << "\t" << bounds[1] << "\t"
                << bounds[2] << "\t" << bounds[3]);
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::GetUnscaledInputBounds(double bounds[4])
{
  bounds[0] = VTK_DOUBLE_MAX;
  bounds[1] = VTK_DOUBLE_MIN;
  bounds[2] = VTK_DOUBLE_MAX;
  bounds[3] = VTK_DOUBLE_MIN;
  for (unsigned int i = 0; i < this->GetNumberOfItems(); i++)
    {
    double b[4];
    this->Internal->Lines[i]->GetUnscaledInputBounds(b);
    if (b[0] < bounds[0]) { bounds[0] = b[0]; }
    if (b[1] > bounds[1]) { bounds[1] = b[1]; }
    if (b[2] < bounds[2]) { bounds[2] = b[2]; }
    if (b[3] > bounds[3]) { bounds[3] = b[3]; }
    }
  vtkDebugMacro(<< "Bounds: " << bounds[0] << "\t" << bounds[1] << "\t"
                << bounds[2] << "\t" << bounds[3]);
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::SetLookupTable(vtkScalarsToColors *lut)
{
  if ( this->LookupTable != lut )
    {
    if ( this->LookupTable)
      {
      this->LookupTable->UnRegister(this);
      }
    this->LookupTable = lut;
    if (lut)
      {
      lut->Register(this);
      }
    this->Modified();
    }
}

//-----------------------------------------------------------------------------
vtkScalarsToColors *vtkPlotFunctionalBag::GetLookupTable()
{
  if ( this->LookupTable == 0 )
    {
    this->CreateDefaultLookupTable();
    }
  return this->LookupTable;
}

//-----------------------------------------------------------------------------
void vtkPlotFunctionalBag::CreateDefaultLookupTable()
{
  if ( this->LookupTable)
    {
    this->LookupTable->UnRegister(this);
    }
  this->LookupTable = vtkLookupTable::New();
  // Consistent Register/UnRegisters.
  this->LookupTable->Register(this);
  this->LookupTable->Delete();
}
