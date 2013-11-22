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

#include "vtkContext2D.h"
#include "vtkPen.h"
#include "vtkBrush.h"
#include "vtkAxis.h"
#include "vtkContextMapper2D.h"
#include "vtkDataArray.h"
#include "vtkTable.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkUnsignedCharArray.h"
#include "vtkLookupTable.h"

#include <algorithm>

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPlotBox);

//-----------------------------------------------------------------------------
vtkPlotBox::vtkPlotBox()
{
  this->BoxWidth = 0.2;
}

//-----------------------------------------------------------------------------
vtkPlotBox::~vtkPlotBox()
{
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

  //painter->ApplyPen(this->Pen);
  this->Brush->SetColor(this->Pen->GetColor());
  vtkNew<vtkPen> blackPen;
  blackPen->SetWidth(this->Pen->GetWidth());
  blackPen->SetColor(0, 0, 0, 128);
  blackPen->SetOpacity(255);    
  painter->ApplyPen(blackPen.GetPointer());
  painter->ApplyBrush(this->Brush);

  // Displays one box plot per row
  vtkTable *table = this->Data->GetInput();

  // Retrieving data array from the input
  vtkDataArray *index = vtkDataArray::SafeDownCast(
    this->Data->GetInputAbstractArrayToProcess(0, this->GetInput()));
  vtkDataArray *qcol[5];
  for (int i = 0; i < 5; i++)
    {
    qcol[i] = vtkDataArray::SafeDownCast(
      this->Data->GetInputAbstractArrayToProcess(i+1, this->GetInput()));
    if (!qcol[i])
      {
      vtkErrorMacro(<< "Input is not correctly set.");
      return false;
      }
    }

  for (int i = 0; i < table->GetNumberOfRows(); ++i)
    {
    // Helper variables for x position
    double x = (!this->UseIndexForXSeries && index) ?
      index->GetVariantValue(i).ToDouble() : i;
    double xpos = x + 0.5 * this->BoxWidth;
    double xneg = x - 0.5 * this->BoxWidth;
    double hBoxW = this->BoxWidth * 0.25;

    double q[5];
    for (int j = 0; j < 5; j++)
      {
      q[j] = qcol[j]->GetVariantValue(i).ToDouble();
      }
    std::sort(q, q+5);

    // Draw the box
    painter->DrawQuad(xpos, q[1], xneg, q[1], xneg, q[3], xpos, q[3]);
    
    // Draw the whiskers
    // Ends of the whiskers match the extremum values of the quartiles
    painter->DrawLine(x, q[0], x, q[1]);
    painter->DrawLine(x, q[3], x, q[4]);
    painter->DrawLine(x - hBoxW, q[0], x + hBoxW, q[0]);
    painter->DrawLine(x - hBoxW, q[4], x + hBoxW, q[4]);

    // Draw the median
    painter->DrawLine(xneg, q[2], xpos, q[2]);
    }
  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::PaintLegend(vtkContext2D *painter, const vtkRectf& rect, int)
{
  painter->ApplyPen(this->Pen);
  painter->DrawLine(rect[0], rect[1] + 0.5 * rect[3],
      rect[0] + rect[2], rect[1] + 0.5 * rect[3]);

  return this->Superclass::PaintLegend(painter, rect, 0);
}

//-----------------------------------------------------------------------------
void vtkPlotBox::GetBounds(double bounds[4])
{
  if (this->Data->GetInput())
    {
    double xRange[2];
    if (!this->UseIndexForXSeries && this->Data->GetInputArrayToProcess(0, this->Data->GetInput()))
      {
      this->Data->
        GetInputArrayToProcess(0, this->Data->GetInput())->GetRange(xRange);
      }
    else
      {
      xRange[0] = 0.;
      xRange[1] = this->Data->GetInputArrayToProcess(1, 
        this->Data->GetInput())->GetNumberOfTuples() - 1;
      }
    bounds[0] = xRange[0] - this->BoxWidth * 0.5;
    bounds[1] = xRange[1] + this->BoxWidth * 0.5;

    if (this->Data->GetInputArrayToProcess(1, this->Data->GetInput()) &&
        this->Data->GetInputArrayToProcess(2, this->Data->GetInput()) &&
        this->Data->GetInputArrayToProcess(3, this->Data->GetInput()) &&
        this->Data->GetInputArrayToProcess(4, this->Data->GetInput()) &&
        this->Data->GetInputArrayToProcess(5, this->Data->GetInput()))
      {
      double yRange[2];
      this->Data->GetInputArrayToProcess(1, this->Data->GetInput())
          ->GetRange(yRange);
      bounds[2] = yRange[0] + this->Pen->GetWidth() * 0.5;
      bounds[3] = yRange[1] + this->Pen->GetWidth() * 0.5;

      for (int i = 1; i < 6; ++i)
        {
        this->Data->GetInputArrayToProcess(i, this->Data->GetInput())
            ->GetRange(yRange);
        if (bounds[2] > yRange[0])
          {
          bounds[2] = yRange[0];
          }
        if (bounds[3] < yRange[1])
          {
          bounds[3] = yRange[1];
          }
        }
      }
    }

  vtkDebugMacro(<< "Bounds: " << bounds[0] << "\t" << bounds[1] << "\t"
                              << bounds[2] << "\t" << bounds[3]);
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetInputData(vtkTable *table)
{
  this->Data->SetInputData(table);
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetInputData(vtkTable *table, const vtkStdString &xColumn,
                              const vtkStdString &yColumn)
{
  // WARN: Do we really need to call the new method with default parameters ?
  // Others classes doesn't... In addition, we assume here that xColumn and
  // yColumn represents the x Position and the q0 values, why ?
  vtkDebugMacro(<< "Setting input, X column = \"" << xColumn.c_str()
                << "\", " << "Y column = \"" << yColumn.c_str() << "\"");

  int sizeToMatch =
      table->GetColumnByName(xColumn.c_str())->GetNumberOfTuples();

  if (table->GetColumn(2)->GetNumberOfTuples() != sizeToMatch ||
      table->GetColumn(3)->GetNumberOfTuples() != sizeToMatch ||
      table->GetColumn(4)->GetNumberOfTuples() != sizeToMatch ||
      table->GetColumn(5)->GetNumberOfTuples() != sizeToMatch )
    {
    vtkErrorMacro(<< "Table is not corretly initialized to use this method");
    return;
    }

  this->SetInputData(table, xColumn, yColumn,
                     table->GetColumnName(2),
                     table->GetColumnName(3),
                     table->GetColumnName(4),
                     table->GetColumnName(5));
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetInputData(vtkTable *table, const vtkStdString &xColumn,
                              const vtkStdString &q0Column,
                              const vtkStdString &q1Column,
                              const vtkStdString &q2Column,
                              const vtkStdString &q3Column,
                              const vtkStdString &q4Column)
{
  vtkDebugMacro(<< "Setting input, X column = \"" << xColumn.c_str()
                << "\", " << "q0 column = \"" << q0Column.c_str() << "\""
                << "\", " << "q1 column = \"" << q1Column.c_str() << "\""
                << "\", " << "q2 column = \"" << q2Column.c_str() << "\""
                << "\", " << "q3 column = \"" << q3Column.c_str() << "\""
                << "\", " << "q4 column = \"" << q4Column.c_str() << "\"");

  this->Data->SetInputData(table);
  this->Data->SetInputArrayToProcess(0, 0, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     xColumn.c_str());
  this->Data->SetInputArrayToProcess(1, 0, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     q0Column.c_str());
  this->Data->SetInputArrayToProcess(2, 0, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     q1Column.c_str());
  this->Data->SetInputArrayToProcess(3, 0, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     q2Column.c_str());
  this->Data->SetInputArrayToProcess(4, 0, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     q3Column.c_str());
  this->Data->SetInputArrayToProcess(5, 0, 0,
                                     vtkDataObject::FIELD_ASSOCIATION_ROWS,
                                     q4Column.c_str());
}

//-----------------------------------------------------------------------------
void vtkPlotBox::SetInputData(vtkTable *table, vtkIdType &xColumn,
                              vtkIdType &q0Column,
                              vtkIdType &q1Column,
                              vtkIdType &q2Column,
                              vtkIdType &q3Column,
                              vtkIdType &q4Column)
{
  this->SetInputData(table,
                     table->GetColumnName(xColumn),
                     table->GetColumnName(q0Column),
                     table->GetColumnName(q1Column),
                     table->GetColumnName(q2Column),
                     table->GetColumnName(q3Column),
                     table->GetColumnName(q4Column));
}

//-----------------------------------------------------------------------------
void vtkPlotBox::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
