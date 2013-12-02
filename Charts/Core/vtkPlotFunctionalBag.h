/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlotFunctionalBag.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkPlotFunctionalBag - Class for drawing an XY line functional 
// bag plot given an input vtkTable.
//
// .SECTION Description
//

#ifndef __vtkPlotFunctionalBag_h
#define __vtkPlotFunctionalBag_h

#include "vtkChartsCoreModule.h" // For export macro
#include "vtkPlot.h"
#include "vtkSmartPointer.h"     // Needed to hold SP ivars
#include "vtkRect.h"             // For vtkRectd ivar

class vtkContext2D;
class vtkContextMapper2D;
class vtkDoubleArray;
class vtkPlotFuntionalBagInternal;
class vtkPoints2D;
class vtkScalarsToColors;
class vtkTable;

class VTKCHARTSCORE_EXPORT vtkPlotFunctionalBag : public vtkPlot
{
public:
  vtkTypeMacro(vtkPlotFunctionalBag, vtkPlot);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Creates a 2D Chart object.
  static vtkPlotFunctionalBag *New();

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
  // This is a convenience function to set the input table for the plot.
  virtual void SetInputData(vtkTable *table);
  
  // Description:
  // This is a convenience function to set the input density table and column for the plot.  
  virtual void SetInputDensityData(vtkTable *table, const vtkStdString &densityColumn, 
    const vtkStdString &variableNameColumn = "");
  virtual void SetInputDensityData(vtkTable *table, vtkIdType densityColumn, 
    vtkIdType variableNameColumn = -1);

  // Description:
  // Specify a lookup table for the mapper to use.
  void SetLookupTable(vtkScalarsToColors *lut);
  vtkScalarsToColors *GetLookupTable();

  // Description:
  // Create default lookup table. Generally used to create one when none
  // is available with the scalar data.
  virtual void CreateDefaultLookupTable();
  
  // Description:
  // Get the plot labels. If this array has a length greater than 1 the index
  // refers to the stacked objects in the plot. See vtkPlotBar for example.
  virtual vtkStringArray *GetLabels();

  // Description:
  // Generate and return the tooltip label string for this plot
  // The segmentIndex parameter is ignored, except for vtkPlotBar
  virtual vtkStdString GetTooltipLabel(const vtkVector2d &plotPos,
                                       vtkIdType seriesIndex,
                                       vtkIdType segmentIndex);

//BTX
  // Description:
  // Function to query a plot for the nearest point to the specified coordinate.
  // Returns the index of the data series with which the point is associated or
  // -1.
  virtual vtkIdType GetNearestPoint(const vtkVector2f& point,
                                    const vtkVector2f& tolerance,
                                    vtkVector2f* location);
//ETX
  // Description:
  // Helper to get the density column array from the density table.
  vtkDoubleArray* GetDensityArray();

protected:
  vtkPlotFunctionalBag();
  ~vtkPlotFunctionalBag();

  // Description:
  // Update the table cache.
  bool UpdateTableCache(vtkTable*);
  
  // Description:
  // Update bag polygons cache.
  bool UpdateBagsCache(vtkTable*, vtkTable*);
  
  // Description:
  // The cache is marked dirty until it has been initialized.
  vtkTimeStamp BuildTime;

  vtkPlotFuntionalBagInternal* Internal;

  // Description:
  // Lookup Table for coloring points by scalar value
  vtkScalarsToColors *LookupTable;

  vtkPoints2D *MedianPoints;
  vtkPoints2D *Q3Points;

  vtkIdType LastNearestSerie;

private:
  vtkPlotFunctionalBag(const vtkPlotFunctionalBag &); // Not implemented.
  void operator=(const vtkPlotFunctionalBag &); // Not implemented.
};

#endif //__vtkPlotFunctionalBag_h
