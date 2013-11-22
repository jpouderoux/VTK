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

// .NAME vtkPlotBox - Class for drawing an a box plot given six columns from
// a vtkTable. Each row of the input table will represent a box. This class
// doesn't compute the quartiles of a data series but directly use them from
// the input tables to display boxes of the vtkPlotBox.
//
// .SECTION Description
//

#ifndef __vtkPlotBox_h
#define __vtkPlotBox_h

#include "vtkChartsCoreModule.h" // For export macro
#include "vtkPlot.h"

class vtkContext2D;
class vtkTable;
class vtkDataArray;
class vtkStdString;


class VTKCHARTSCORE_EXPORT vtkPlotBox : public vtkPlot
{
public:
  vtkTypeMacro(vtkPlotBox, vtkPlot);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Creates a 2D Chart object.
  static vtkPlotBox *New();

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
  // Get the bounds for this box plot as (Xmin, Xmax, Ymin, Ymax).
  virtual void GetBounds(double bounds[4]);

  // Description:
  // Set the input, we are expecting a vtkTable with six columns. The first
  // column represents the x position of a displayed box. The five others
  // columns represent the quartiles used to display the box.
  // Inherited method will call the last SetInputData method with default
  // paramaters.
  virtual void SetInputData(vtkTable *table);
  virtual void SetInputData(vtkTable *table, const vtkStdString &xColumn,
                            const vtkStdString &yColumn);
  virtual void SetInputData(vtkTable *table, const vtkStdString &xColumn,
                            const vtkStdString &q0Column,
                            const vtkStdString &q1Column,
                            const vtkStdString &q2Column,
                            const vtkStdString &q3Column,
                            const vtkStdString &q4Column);

  void SetInputData(vtkTable *table, vtkIdType &xColumn,
                    vtkIdType &q0Column,
                    vtkIdType &q1Column,
                    vtkIdType &q2Column,
                    vtkIdType &q3Column,
                    vtkIdType &q4Column);


  // Description:
  // Set the width of the box.
  vtkSetMacro(BoxWidth, float);

  // Description:
  // Get the width of the box.
  vtkGetMacro(BoxWidth, float);

  //BTX
protected:
  vtkPlotBox();
  ~vtkPlotBox();

  float BoxWidth;

  // Description:
  // The point cache is marked dirty until it has been initialized.
  vtkTimeStamp BuildTime;

private:
  vtkPlotBox(const vtkPlotBox &); // Not implemented.
  void operator=(const vtkPlotBox &); // Not implemented.

  //ETX
};

#endif //__vtkPlotBox_h
