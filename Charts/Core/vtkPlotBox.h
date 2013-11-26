/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlotBox.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkPlotBox - Class for drawing a box plot given one column of
// a vtkTable.
//
// .SECTION Description
// The selected column of the table must contains the quartiles and the
// extremum values in the 5 first lines (the order does not matter).
//
// .SECTION See also
// vtlPlotPoints vtkComputeQuartiles
//
// .SECTION Thanks
// This class was written by Kitware SAS and supported by EDF - www.edf.fr

#ifndef __vtkPlotBox_h
#define __vtkPlotBox_h

#include "vtkChartsCoreModule.h" // For export macro
#include "vtkPlotPoints.h"

class VTKCHARTSCORE_EXPORT vtkPlotBox : public vtkPlotPoints
{
public:
  vtkTypeMacro(vtkPlotBox, vtkPlotPoints);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Creates a 2D Chart object.
  static vtkPlotBox *New();

  // Description:
  // Perform any updates to the item that may be necessary before rendering.
  // The scene should take care of calling this on all items before their
  // Paint function is invoked.
  virtual void Update();

  // Description:
  // Paint event for the XY plot, called whenever the chart needs to be drawn.
  virtual bool Paint(vtkContext2D *painter);

  // Description:
  // Paint legend event for the XY plot, called whenever the legend needs the
  // plot items symbol/mark/line drawn. A rect is supplied with the lower left
  // corner of the rect (elements 0 and 1) and with width x height (elements 2
  // and 3). The plot can choose how to fill the space supplied.
  virtual bool PaintLegend(vtkContext2D *painter, const vtkRectf& rect,
                           int legendIndex);

  // Description:
  // Get the bounds for this plot as (Xmin, Xmax, Ymin, Ymax).
  virtual void GetBounds(double bounds[4]);

  // Description:
  // Get the non-log-scaled bounds on chart inputs for this plot as (Xmin, Xmax, Ymin, Ymax).
  virtual void GetUnscaledInputBounds(double bounds[4]);

  // Description:
  // Get/Set the width of the box. Default is 0.2
  vtkGetMacro(BoxWidth, float);
  vtkSetMacro(BoxWidth, float);

//BTX
protected:
  vtkPlotBox();
  ~vtkPlotBox();

  float BoxWidth;
  int Position;

private:
  vtkPlotBox(const vtkPlotBox &); // Not implemented.
  void operator=(const vtkPlotBox &); // Not implemented.
//ETX
};

#endif //__vtkPlotLine_h
