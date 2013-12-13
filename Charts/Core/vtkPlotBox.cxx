/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlotBox.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPlotBox.h"

#include "vtkChartBox.h"
#include "vtkContext2D.h"
#include "vtkAxis.h"
#include "vtkBrush.h"
#include "vtkPen.h"
#include "vtkNew.h"
#include "vtkDoubleArray.h"
#include "vtkVector.h"
#include "vtkContextDevice2D.h"
#include "vtkContextMapper2D.h"
#include "vtkTable.h"
#include "vtkDataArray.h"
#include "vtkIdTypeArray.h"
#include "vtkTimeStamp.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkUnsignedCharArray.h"
#include "vtkStringArray.h"

// Need to turn some arrays of strings into categories
#include "vtkStringToCategory.h"

#include "vtkObjectFactory.h"

#include <vector>
#include <algorithm>

class vtkPlotBox::Private :
    public std::vector< std::vector<double> >
{
public:
  Private()
  {
    this->SelectionInitialized = false;
  }

  bool SelectionInitialized;
};


//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPlotBox)

//-----------------------------------------------------------------------------
vtkPlotBox::vtkPlotBox()
{
  this->Storage = new vtkPlotBox::Private;
  this->Pen->SetColor(255, 0, 0, 25);
  this->BoxWidth = 20.;
  this->LookupTable = 0;
  this->Colors = 0;
  this->ScalarVisibility = 0;
}

//-----------------------------------------------------------------------------
vtkPlotBox::~vtkPlotBox()
{
  delete this->Storage;

  if (this->LookupTable)
    {
    this->LookupTable->UnRegister(this);
    }
  if (this->Colors)
    {
    this->Colors->UnRegister(this);
    }
}

//-----------------------------------------------------------------------------
void vtkPlotBox::Update()
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

  this->UpdateTableCache(table);
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::Paint(vtkContext2D *painter)
{
  // This is where everything should be drawn, or dispatched to other methods.
  vtkDebugMacro(<< "Paint event called in vtkPlotBox.");

  if (!this->Visible)
    {
    return false;
    }

  painter->ApplyPen(this->Pen);

  if (this->Storage->size() == 0)
    {
    return false;
    }

  size_t cols = this->Storage->size();
  size_t rows = this->Storage->at(0).size();

  vtkIdType selection = 0;
  vtkIdType id = 0;
  vtkIdType selectionSize = 0;
  if (this->Selection)
    {
    selectionSize = this->Selection->GetNumberOfTuples();
    if (selectionSize)
      {
      this->Selection->GetTupleValue(selection, &id);
      }
    }

    // Update the axis positions
  vtkChartBox *parent = vtkChartBox::SafeDownCast(this->Parent);

  if (cols > 1)
    {
    this->BoxWidth = (parent->GetAxis(1)->GetPoint1()[0] -
     parent->GetAxis(0)->GetPoint1()[0]) / 2.;
    }

  for (int i = 0; i < cols; i++)
    {
    vtkNew<vtkBrush> brush;
    brush->SetColor(this->Pen->GetColor());
    vtkNew<vtkPen> blackPen;
    blackPen->SetWidth(this->Pen->GetWidth());
    blackPen->SetColor(0, 0, 0, 128);
    blackPen->SetOpacity(255);
    painter->ApplyPen(blackPen.GetPointer());
    painter->ApplyBrush(brush.GetPointer());

    // Helper variables for x position
    double x = parent->GetAxis(int(i))->GetPoint1()[0];
    double xpos = x + 0.5 * this->BoxWidth;
    double xneg = x - 0.5 * this->BoxWidth;
    double hBoxW = this->BoxWidth * 0.25;

    // Fetch the quartiles and median
    double q[5];
    for (int j = 0; j < 5; j++)
      {
      q[j] = (*this->Storage)[i][j];
      }
    std::sort(q, q+5);

    // Draw the box
    painter->DrawQuad(xpos, q[1], xneg, q[1], xneg, q[3], xpos, q[3]);

    // Draw the whiskers: ends of the whiskers match the
    // extremum values of the quartiles
    painter->DrawLine(x, q[0], x, q[1]);
    painter->DrawLine(x - hBoxW, q[0], x + hBoxW, q[0]);
    painter->DrawLine(x, q[3], x, q[4]);
    painter->DrawLine(x - hBoxW, q[4], x + hBoxW, q[4]);

    // Draw the median
    vtkNew<vtkPen> whitePen;
    unsigned char brushColor[4];
    brush->GetColor(brushColor);
    // Use a gray pen if the brush is black so the median is always visible
    if (brushColor[0] == 0 && brushColor[1] == 0 && brushColor[2] == 0)
      {
      whitePen->SetWidth(std::max(1.0, this->Pen->GetWidth() - 1.0));
      whitePen->SetColor(128, 128, 128, 128);
      whitePen->SetOpacity(this->Pen->GetOpacity());
      painter->ApplyPen(whitePen.GetPointer());
      }
    painter->DrawLine(xneg, q[2], xpos, q[2]);
  }

  // Draw all of the lines
  /*painter->ApplyPen(this->Pen);
  int ncComps(0);
  if (this->ScalarVisibility && this->Colors)
    {
    ncComps = static_cast<int>(this->Colors->GetNumberOfComponents());
    }
  if (this->ScalarVisibility && this->Colors && ncComps == 4)
    {
    for (size_t i = 0, nc = 0; i < rows; ++i, nc += ncComps)
      {
      for (size_t j = 0; j < cols; ++j)
        {
        line[j].Set(this->Storage->AxisPos[j], (*this->Storage)[j][i]);
        }
      painter->GetPen()->SetColor(this->Colors->GetPointer(nc));
      painter->DrawPoly(line[0].GetData(), static_cast<int>(cols));
      }
    }
  else
    {
    for (size_t i = 0; i < rows; ++i)
      {
      for (size_t j = 0; j < cols; ++j)
        {
        line[j].Set(this->Storage->AxisPos[j], (*this->Storage)[j][i]);
        }
      painter->DrawPoly(line[0].GetData(), static_cast<int>(cols));
      }
    }

  // Now draw the selected lines
  if (this->Selection)
    {
    painter->GetPen()->SetColor(255, 0, 0, 100);
    for (vtkIdType i = 0; i < this->Selection->GetNumberOfTuples(); ++i)
      {
      for (size_t j = 0; j < cols; ++j)
        {
        this->Selection->GetTupleValue(i, &id);
        line[j].Set(this->Storage->AxisPos[j], (*this->Storage)[j][id]);
        }
      painter->DrawPoly(line[0].GetData(), static_cast<int>(cols));
      }
    }

  delete[] line;*/

  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::PaintLegend(vtkContext2D *painter,
                                             const vtkRectf& rect, int)
{
  painter->ApplyPen(this->Pen);
  painter->DrawLine(rect[0]          , rect[1] + 0.5 * rect[3],
                    rect[0] + rect[2], rect[1] + 0.5 * rect[3]);
  return true;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::GetBounds(double *)
{

}

//-----------------------------------------------------------------------------
bool vtkPlotBox::SetSelectionRange(int axis, float low,
                                                   float high)
{
  if (!this->Selection)
    {
    return false;
    }
  if (this->Storage->SelectionInitialized)
    {
    // Further refine the selection that has already been made
    vtkIdTypeArray *array = vtkIdTypeArray::New();
    std::vector<double>& col = this->Storage->at(axis);
    for (vtkIdType i = 0; i < this->Selection->GetNumberOfTuples(); ++i)
      {
      vtkIdType id = 0;
      this->Selection->GetTupleValue(i, &id);
      if (col[id] >= low && col[id] <= high)
        {
        // Remove this point - no longer selected
        array->InsertNextValue(id);
        }
      }
    this->Selection->DeepCopy(array);
    array->Delete();
    }
  else
    {
    // First run - ensure the selection list is empty and build it up
    std::vector<double>& col = this->Storage->at(axis);
    for (size_t i = 0; i < col.size(); ++i)
      {
      if (col[i] >= low && col[i] <= high)
        {
        // Remove this point - no longer selected
        this->Selection->InsertNextValue(i);
        }
      }
    this->Storage->SelectionInitialized = true;
    }
  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::ResetSelectionRange()
{
  this->Storage->SelectionInitialized = false;
  if (this->Selection)
    {
    this->Selection->SetNumberOfTuples(0);
    }
  return true;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetInputData(vtkTable* table)
{
  if (table == this->Data->GetInput() && (!table ||
                                          table->GetMTime() < this->BuildTime))
    {
    return;
    }

  bool updateVisibility = table != this->Data->GetInput();
  this->vtkPlot::SetInputData(table);
  vtkChartBox *parent =
      vtkChartBox::SafeDownCast(this->Parent);

  if (parent && table && updateVisibility)
    {
    parent->SetColumnVisibilityAll(false);
    // By default make the first 10 columns visible in a plot.
    for (vtkIdType i = 0; i < table->GetNumberOfColumns() && i < 10; ++i)
      {
      parent->SetColumnVisibility(table->GetColumnName(i), true);
      }
    }
  else if (parent && updateVisibility)
    {
    // No table, therefore no visible columns
    parent->GetVisibleColumns()->SetNumberOfTuples(0);
    }
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::UpdateTableCache(vtkTable *table)
{
  // Each axis is a column in our storage array, they are scaled from 0.0 to 1.0
  vtkChartBox *parent =
      vtkChartBox::SafeDownCast(this->Parent);
  if (!parent || !table || table->GetNumberOfColumns() == 0)
    {
    return false;
    }

  vtkStringArray* cols = parent->GetVisibleColumns();

  this->Storage->resize(cols->GetNumberOfTuples());
  vtkIdType rows = table->GetNumberOfRows();

  for (vtkIdType i = 0; i < cols->GetNumberOfTuples(); ++i)
    {
    std::vector<double>& col = this->Storage->at(i);
    col.resize(rows);
    vtkSmartPointer<vtkDataArray> data =
        vtkDataArray::SafeDownCast(table->GetColumnByName(cols->GetValue(i)));
    if (!data)
      {
      continue;
      }

    vtkAxis* axis = parent->GetAxis(i);
    // Also need the range from the appropriate axis, to normalize points
    double min = axis->GetUnscaledMinimum();
    double max = axis->GetUnscaledMaximum();
    double scale = 1.0f / (max - min);

    for (vtkIdType j = 0; j < rows; ++j)
      {
      col[j] = (data->GetTuple1(j) - min) * scale;
      }
    }

  // Additions for color mapping
  if (this->ScalarVisibility && !this->ColorArrayName.empty())
    {
    vtkDataArray* c =
      vtkDataArray::SafeDownCast(table->GetColumnByName(this->ColorArrayName));
    // TODO: Should add support for categorical coloring & try enum lookup
    if (c)
      {
      if (!this->LookupTable)
        {
        this->CreateDefaultLookupTable();
        }
      this->Colors = this->LookupTable->MapScalars(c, VTK_COLOR_MODE_MAP_SCALARS, -1);
      // Consistent register and unregisters
      this->Colors->Register(this);
      this->Colors->Delete();
      }
    else
      {
      this->Colors->UnRegister(this);
      this->Colors = 0;
      }
    }

  this->BuildTime.Modified();
  return true;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetLookupTable(vtkScalarsToColors *lut)
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
vtkScalarsToColors *vtkPlotBox::GetLookupTable()
{
  if ( this->LookupTable == 0 )
    {
    this->CreateDefaultLookupTable();
    }
  return this->LookupTable;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::CreateDefaultLookupTable()
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

//-----------------------------------------------------------------------------
void vtkPlotBox::SelectColorArray(const vtkStdString &arrayName)
{
  vtkTable *table = this->Data->GetInput();
  if (!table)
    {
    vtkDebugMacro(<< "SelectColorArray called with no input table set.");
    return;
    }
  if (this->ColorArrayName == arrayName)
    {
    return;
    }
  for (vtkIdType c = 0; c < table->GetNumberOfColumns(); ++c)
    {
    if (table->GetColumnName(c) == arrayName)
      {
      this->ColorArrayName = arrayName;
      this->Modified();
      return;
      }
    }
  vtkDebugMacro(<< "SelectColorArray called with invalid column name.");
  this->ColorArrayName = "";
  this->Modified();
  return;
}

//-----------------------------------------------------------------------------
vtkStdString vtkPlotBox::GetColorArrayName()
{
  return this->ColorArrayName;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SelectColorArray(vtkIdType arrayNum)
{
  vtkTable *table = this->Data->GetInput();
  if (!table)
    {
    vtkDebugMacro(<< "SelectColorArray called with no input table set.");
    return;
    }
  vtkDataArray *col = vtkDataArray::SafeDownCast(table->GetColumn(arrayNum));
  // TODO: Should add support for categorical coloring & try enum lookup
  if (!col)
    {
    vtkDebugMacro(<< "SelectColorArray called with invalid column index");
    return;
    }
  else
    {
    if (this->ColorArrayName == table->GetColumnName(arrayNum))
      {
      return;
      }
    else
      {
      this->ColorArrayName = table->GetColumnName(arrayNum);
      this->Modified();
      }
    }
}

//-----------------------------------------------------------------------------
void vtkPlotBox::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
