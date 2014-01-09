/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRGrid.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkRGrid - topologically regular array of data
// .SECTION Description
// vtkRGrid is a data object that is a concrete implementation of
// vtkDataSet. vtkRGrid represents a geometric structure that is a
// topologically regular array of points. The topology is that of a cube that
// has been subdivided into a regular array of smaller cubes. Each point/cell
// can be addressed with i-j-k indices. Examples include finite difference
// grids.
//
// The order and number of points must match that specified by the dimensions
// of the grid. The point order increases in i fastest (from 0<=i<dims[0]),
// then j (0<=j<dims[1]), then k (0<=k<dims[2]) where dims[] are the
// dimensions of the grid in the i-j-k topological directions. The number of
// points is dims[0]*dims[1]*dims[2]. The same is true for the cells of the
// grid. The order and number of cells must match that specified by the
// dimensions of the grid. The cell order increases in i fastest (from
// 0<=i<(dims[0]-1)), then j (0<=j<(dims[1]-1)), then k (0<=k<(dims[2]-1))
// The number of cells is (dims[0]-1)*(dims[1]-1)*(dims[2]-1).
//
// A unusual feature of vtkRGrid is the ability to blank,
// or "turn-off" points and cells in the dataset. This is controlled by
// defining a "blanking array" whose values (0,1) specify whether
// a point should be blanked or not.

#ifndef __vtkRGrid_h
#define __vtkRGrid_h

#include "vtkCommonDataModelModule.h" // For export macro
#include "vtkPointSet.h"

#include "vtkRGrid.h" // Needed for inline methods

class vtkEmptyCell;
class vtkHexahedron;
class vtkUnsignedCharArray;

class VTKCOMMONDATAMODEL_EXPORT vtkRGrid : public vtkPointSet
{
public:
  static vtkRGrid *New();

  vtkTypeMacro(vtkRGrid,vtkPointSet);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Return what type of dataset this is.
  int GetDataObjectType() {return VTK_R_GRID;}

  // Description:
  // Copy the geometric and topological structure of an input poly data object.
  void CopyStructure(vtkDataSet *ds);

  // Description:
  // Standard vtkDataSet API methods. See vtkDataSet for more information.
  vtkIdType GetNumberOfPoints() {return vtkPointSet::GetNumberOfPoints();}
  double *GetPoint(vtkIdType ptId) {return this->vtkPointSet::GetPoint(ptId);}
  void GetPoint(vtkIdType ptId, double p[3])
    {this->vtkPointSet::GetPoint(ptId,p);}
  vtkCell *GetCell(vtkIdType cellId);
  void GetCell(vtkIdType cellId, vtkGenericCell *cell);
  void GetCellBounds(vtkIdType cellId, double bounds[6]);
  int GetCellType(vtkIdType cellId);
  vtkIdType GetNumberOfCells();
  void GetCellPoints(vtkIdType cellId, vtkIdList *ptIds);
  void GetPointCells(vtkIdType ptId, vtkIdList *cellIds);

  void Initialize();
  int GetMaxCellSize() {return 8;}; //hexahedron is the largest
  void GetCellNeighbors(vtkIdType cellId, vtkIdList *ptIds,
                        vtkIdList *cellIds);

  // Description:
  // following methods are specific to structured grid
  void SetDimensions(int i, int j, int k);
  void SetDimensions(int dim[3]);

  // Description:
  // Get dimensions of this structured points dataset.
  virtual int *GetDimensions ();
  virtual void GetDimensions (int dim[3]);

  // Description:
  // Return the dimensionality of the data.
  int GetDataDimension() { return 3; }

  // Description:
  // Different ways to set the extent of the data array.  The extent
  // should be set before the "Scalars" are set or allocated.
  // The Extent is stored  in the order (X, Y, Z).
  void SetExtent(int extent[6]);
  void SetExtent(int x1, int x2, int y1, int y2, int z1, int z2);
  vtkGetVector6Macro(Extent, int);

  // Description:
  // Return the actual size of the data in kilobytes. This number
  // is valid only after the pipeline has updated. The memory size
  // returned is guaranteed to be greater than or equal to the
  // memory required to represent the data (e.g., extra space in
  // arrays, etc. are not included in the return value). THIS METHOD
  // IS THREAD SAFE.
  unsigned long GetActualMemorySize();

  // Description:
  // Shallow and Deep copy.
  void ShallowCopy(vtkDataObject *src);
  void DeepCopy(vtkDataObject *src);

  // Description:
  // The extent type is a 3D extent
  int GetExtentType() { return VTK_3D_EXTENT; }

  // Description:
  // Given the node dimensions of this grid instance, this method computes the
  // node dimensions. The value in each dimension can will have a lowest value
  // of "1" such that computing the total number of cells can be achieved by
  // simply by cellDims[0]*cellDims[1]*cellDims[2].
  void GetCellDims( int cellDims[3] );

  // Description:
  // Reallocates and copies to set the Extent to the UpdateExtent.
  // This is used internally when the exact extent is requested,
  // and the source generated more than the update extent.
  virtual void Crop(const int* updateExtent);

  //BTX
  // Description:
  // Retrieve an instance of this class from an information object.
  static vtkRGrid* GetData(vtkInformation* info);
  static vtkRGrid* GetData(vtkInformationVector* v, int i=0);
  //ETX

  // Description:
  // Get a point in the grid. If adjustForExtent is true, (i,j,k) is
  // interpreted as a position relative to the beginning of the extent.
  // If adjustForExtent is false, (i,j,k) is interpreted literally
  // and the (i,j,k) point of the grid is returned regardless of the
  // extent beginning.
  // The point coordinate is returned in 'p'.
  // The default adjustForExtent is true.
  void GetPoint(int i, int j, int k, double p[3], bool adjustForExtent = true);

protected:
  vtkRGrid();
  ~vtkRGrid();

  // for the GetCell method
  vtkHexahedron *Hexahedron;
  vtkEmptyCell *EmptyCell;

  int Dimensions[3];

  int Extent[6];

  // Description:
  // Compute the range of the scalars and cache it into ScalarRange
  // only if the cache became invalid (ScalarRangeComputeTime).
  virtual void ComputeScalarRange();

private:
  // Description:
  // For legacy compatibility. Do not use.
  void GetCellNeighbors(vtkIdType cellId, vtkIdList& ptIds, vtkIdList& cellIds)
    {this->GetCellNeighbors(cellId, &ptIds, &cellIds);}

  // Internal method used by DeepCopy and ShallowCopy.
  void InternalRGridCopy(vtkRGrid *src);

private:
  vtkRGrid(const vtkRGrid&);  // Not implemented.
  void operator=(const vtkRGrid&);  // Not implemented.
};


inline vtkIdType vtkRGrid::GetNumberOfCells()
{
  int nCells=1;
  int dims[3];
  int i;

  this->GetDimensions(dims);
  for (i=0; i<3; i++)
    {
    if (dims[i] <= 0)
      {
      return 0;
      }
    if (dims[i] > 1)
      {
      nCells *= (dims[i]-1);
      }
    }

  return nCells;
}

#endif
