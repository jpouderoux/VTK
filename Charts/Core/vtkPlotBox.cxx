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

#include "vtkBrush.h"
#include "vtkContext2D.h"
#include "vtkObjectFactory.h"
#include "vtkPen.h"
#include "vtkRect.h"
#include "vtkPoints2D.h"
#include "vtkTable.h"

#include <algorithm>

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPlotBox);

//-----------------------------------------------------------------------------
vtkPlotBox::vtkPlotBox()
{
  this->MarkerStyle = vtkPlotPoints::NONE;
  this->Position = 0;
  this->BoxWidth = 0.2;
  this->TooltipDefaultLabelFormat = "%l: %y";
}

//-----------------------------------------------------------------------------
vtkPlotBox::~vtkPlotBox()
{
}

//-----------------------------------------------------------------------------
void vtkPlotBox::Update()
{
  this->Superclass::Update();

  // Compute the x-position of the plot by fetching the plot index in the
  // parent's visible children.
  if (this->Points && this->Points->GetNumberOfPoints() >= 5 && this->Visible)
    {
    vtkAbstractContextItem *parent = this->GetParent();
    unsigned int nbItems = parent->GetNumberOfItems();
    unsigned int itemPos = 0;
    for (unsigned int i = 0; i < nbItems; i++)
      {
      if (parent->GetItem(i) == this)
        {
        this->Position = itemPos;
        // Set the points x-position for other plot functionalities
        // like selection, etc.
        for (int j = 0; j < 5; j++)
         {
           this->Points->SetPoint(j,
             this->Position, this->Points->GetPoint(j)[1]);
         }
        break;
        }
      else
        {
        itemPos += parent->GetItem(i)->GetVisible() ? 1 : 0;
        }
      }
    }
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::Paint(vtkContext2D *painter)
{
  // This is where everything should be drawn, or dispatched to other methods.
  vtkDebugMacro(<< "Paint event called in vtkPlotBox.");

  if (!this->Visible || this->Points->GetNumberOfPoints() < 5 || !this->Points)
    {
    return false;
    }

  this->Brush->SetColor(this->Pen->GetColor());
  vtkNew<vtkPen> blackPen;
  blackPen->SetWidth(this->Pen->GetWidth());
  blackPen->SetColor(0, 0, 0, 128);
  blackPen->SetOpacity(255);
  painter->ApplyPen(blackPen.GetPointer());
  painter->ApplyBrush(this->Brush);

  // Helper variables for x position
  double x = this->Position;
  double xpos = x + 0.5 * this->BoxWidth;
  double xneg = x - 0.5 * this->BoxWidth;
  double hBoxW = this->BoxWidth * 0.25;

  // Fetch the quartiles and median
  double q[5];
  for (int j = 0; j < 5; j++)
    {
    q[j] = this->Points->GetPoint(j)[1];
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
  this->Brush->GetColor(brushColor);
  // Use a gray pen if the brush is black so the median is always visible
  if (brushColor[0] == 0 && brushColor[1] == 0 && brushColor[2] == 0)
    {
    whitePen->SetWidth(std::max(1.0, this->Pen->GetWidth() - 1.0));
    whitePen->SetColor(128, 128, 128, 128);
    whitePen->SetOpacity(this->Pen->GetOpacity());
    painter->ApplyPen(whitePen.GetPointer());
    }
  painter->DrawLine(xneg, q[2], xpos, q[2]);

  return true;
}

//-----------------------------------------------------------------------------
bool vtkPlotBox::PaintLegend(vtkContext2D *painter, const vtkRectf& rect, int)
{
  vtkNew<vtkPen> blackPen;
  blackPen->SetWidth(1.0);
  blackPen->SetColor(0, 0, 0, 128);
  blackPen->SetOpacity(255);
  painter->ApplyPen(blackPen.GetPointer());
  painter->ApplyBrush(this->Brush);
  painter->DrawRect(rect[0], rect[1], rect[2], rect[3]);

  return this->Superclass::PaintLegend(painter, rect, 0);
}

//-----------------------------------------------------------------------------
void vtkPlotBox::GetBounds(double bounds[4])
{
  // We can use the BadPoints array to skip the bad points
  if (!this->Points)
    {
    return;
    }

  this->GetUnscaledInputBounds(bounds);

  if (this->LogX)
    {
    bounds[0] = log10(bounds[0]);
    bounds[1] = log10(bounds[1]);
    }
  if (this->LogY)
    {
    bounds[2] = log10(bounds[2]);
    bounds[3] = log10(bounds[3]);
    }

  vtkDebugMacro(
    << "Bounds: " << bounds[0] << "\t" << bounds[1] << "\t"
    << bounds[2] << "\t" << bounds[3]);
}

//-----------------------------------------------------------------------------
void vtkPlotBox::GetUnscaledInputBounds(double bounds[4])
{
  // We can use the BadPoints array to skip the bad points
  if (!this->Points)
    {
    return;
    }

  // Get q0 & q4 (min and max)
  double q0 = this->Points->GetPoint(0)[1];
  double q4 = q0;
  for (int j = 1; j < 5; j++)
    {
    double q = this->Points->GetPoint(j)[1];
    if (q < q0)
      {
      q0 = q;
      }
    else if (q > q4)
      {
      q4 = q;
      }
    }

  bounds[0] = this->Position - 0.5 * this->BoxWidth;
  bounds[1] = this->Position + 0.5 * this->BoxWidth;
  bounds[2] = q0;
  bounds[3] = q4;

  vtkDebugMacro(
    << "Bounds: " << bounds[0] << "\t" << bounds[1] << "\t"
    << bounds[2] << "\t" << bounds[3]);
}
//-----------------------------------------------------------------------------
void vtkPlotBox::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
