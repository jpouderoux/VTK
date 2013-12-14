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
#include "vtkDataSetAttributes.h"
#include "vtkAbstractArray.h"

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
  this->ScalarVisibility = 0;
  this->TooltipDefaultLabelFormat = "%y";
}

//-----------------------------------------------------------------------------
vtkPlotBox::~vtkPlotBox()
{
  delete this->Storage;

  if (this->LookupTable)
    {
    this->LookupTable->UnRegister(this);
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

  if (this->Storage->size() == 0 || this->Storage->at(0).size() < 5)
    {
    return false;
    }

  size_t cols = this->Storage->size();

  vtkChartBox *parent = vtkChartBox::SafeDownCast(this->Parent);

  for (int i = 0; i < cols; i++)
    {
    vtkStdString colName = parent->GetVisibleColumns()->GetValue(i);
    int index;
    this->GetInput()->GetRowData()->GetAbstractArray(colName.c_str(), index);
    double rgb[4];
    this->LookupTable->GetIndexedColor(index, rgb);
    vtkNew<vtkBrush> brush;
    brush->SetColor(rgb[0] * 255., rgb[1] * 255, rgb[2] * 255, 255);
    vtkNew<vtkPen> blackPen;
    blackPen->SetWidth(this->Pen->GetWidth());
    blackPen->SetColor(0, 0, 0, 128);
    blackPen->SetOpacity(255);

    if (parent->GetSelectedColumn() == i)
      {
      unsigned char* col = brush->GetColor();
      brush->SetColor(col[0]^255, col[1]^255, col[2]^255, 255);
      }

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

  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::PaintLegend(vtkContext2D *painter,
                             const vtkRectf& rect, int)
{
  return true;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::GetBounds(double *)
{
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetInputData(vtkTable* table)
{
  if (table == this->Data->GetInput() &&
    (!table || table->GetMTime() < this->BuildTime))
    {
    return;
    }

  bool updateVisibility = table != this->Data->GetInput();
  this->vtkPlot::SetInputData(table);
  vtkChartBox *parent = vtkChartBox::SafeDownCast(this->Parent);

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
  if (!this->LookupTable)
    {
    this->CreateDefaultLookupTable();
    }
}

namespace
{
// See if the point is within tolerance.
bool inRange(const vtkVector2f& point, const vtkVector2f& tol,
             const vtkVector2f& current)
{
  return current.GetX() > point.GetX() - tol.GetX() &&
         current.GetX() < point.GetX() + tol.GetX() &&
         current.GetY() > point.GetY() - tol.GetY() &&
         current.GetY() < point.GetY() + tol.GetY();
}
}

//-----------------------------------------------------------------------------
vtkIdType vtkPlotBox::GetNearestPoint(const vtkVector2f& point,
                                      const vtkVector2f& tol,
                                      vtkVector2f* location)
{
  size_t cols = this->Storage->size();

  vtkChartBox *parent = vtkChartBox::SafeDownCast(this->Parent);

  for (int i = 0; i < cols; i++)
    {
    vtkVector2f v;
    v.SetX(parent->GetAxis(int(i))->GetPoint1()[0]);
    for (int j = 0; j < 5; j++)
      {
      v.SetY((*this->Storage)[i][j]);
      if (inRange(point, tol, v))
        {
        vtkAxis* axis = parent->GetAxis(i);
        double min = axis->GetUnscaledMinimum();
        double max = axis->GetUnscaledMaximum();
        double scale = 1.0f / (max - min);
        double y = (*this->Storage)[i][j] / scale + min;
        location->SetX(i);
        location->SetY(y);
        return static_cast<int>(i);
        }
      }
    }
  return -1;

}
//-----------------------------------------------------------------------------
bool vtkPlotBox::UpdateTableCache(vtkTable *table)
{
  // Each boxplot is a column in our storage array,
  // they are scaled from 0.0 to 1.0
  vtkChartBox *parent = vtkChartBox::SafeDownCast(this->Parent);

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

  this->BuildTime.Modified();
  return true;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetLookupTable(vtkScalarsToColors *lut)
{
  if (this->LookupTable != lut)
    {
    if (this->LookupTable)
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
  if (this->LookupTable == 0)
    {
    this->CreateDefaultLookupTable();
    }
  return this->LookupTable;
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetColumnColor(const vtkStdString& colName, double *rgb)
{
  if (this->LookupTable == 0)
    {
    this->CreateDefaultLookupTable();
    }
  int index;
  this->GetInput()->GetRowData()->GetAbstractArray(colName.c_str(), index);
  vtkLookupTable* lut = vtkLookupTable::SafeDownCast(this->LookupTable);
  if (index >= 0 && lut)
    {
    lut->SetTableValue(index, rgb[0], rgb[1], rgb[2]);
    lut->Build();
    }
}

//-----------------------------------------------------------------------------
void vtkPlotBox::CreateDefaultLookupTable()
{
  if (this->LookupTable)
    {
    this->LookupTable->UnRegister(this);
    }
  vtkLookupTable* lut = vtkLookupTable::New();
  this->LookupTable = lut;
  // Consistent Register/UnRegisters.
  this->LookupTable->Register(this);
  this->LookupTable->Delete();
  vtkTable *table = this->GetInput();
  lut->SetNumberOfColors(table->GetNumberOfColumns());
  this->LookupTable->Build();
}

//-----------------------------------------------------------------------------
void vtkPlotBox::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
