/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRGrid.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkRGrid.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLinks.h"
#include "vtkEmptyCell.h"
#include "vtkGenericCell.h"
#include "vtkHexahedron.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkQuad.h"
#include "vtkStructuredData.h"
#include "vtkStructuredVisibilityConstraint.h"
#include "vtkVertex.h"

vtkStandardNewMacro(vtkRGrid);

vtkSetObjectImplementationMacro(vtkRGrid, FacesConnectivityMaskArray, vtkUnsignedCharArray);

#define vtkAdjustBoundsMacro(A, B) \
  A[0] = (B[0] < A[0] ? B[0] : A[0]);   A[1] = (B[0] > A[1] ? B[0] : A[1]); \
  A[2] = (B[1] < A[2] ? B[1] : A[2]);   A[3] = (B[1] > A[3] ? B[1] : A[3]); \
  A[4] = (B[2] < A[4] ? B[2] : A[4]);   A[5] = (B[2] > A[5] ? B[2] : A[5])

//----------------------------------------------------------------------------
vtkRGrid::vtkRGrid()
{
  this->Hexahedron = vtkHexahedron::New();
  this->EmptyCell = vtkEmptyCell::New();

  this->Cells = 0;
  this->Links = 0;

  this->Dimensions[0] = 0;
  this->Dimensions[1] = 0;
  this->Dimensions[2] = 0;

  this->FacesConnectivityMaskArray = 0;
  this->CellVisibility = vtkStructuredVisibilityConstraint::New();

  int extent[6] = { 0, -1, 0, -1, 0, -1 };
  memcpy(this->Extent, extent, 6 * sizeof(int));

  this->Information->Set(vtkDataObject::DATA_EXTENT_TYPE(), VTK_3D_EXTENT);
  this->Information->Set(vtkDataObject::DATA_EXTENT(), this->Extent, 6);
}

//----------------------------------------------------------------------------
vtkRGrid::~vtkRGrid()
{
  this->Hexahedron->Delete();
  this->EmptyCell->Delete();
  this->CellVisibility->Delete();
  if (this->FacesConnectivityMaskArray)
    {
    this->FacesConnectivityMaskArray->UnRegister(this);
    }
}

//----------------------------------------------------------------------------
// Copy the geometric and topological structure of an input structured grid.
void vtkRGrid::CopyStructure(vtkDataSet *ds)
{
  vtkRGrid *sg = static_cast<vtkRGrid*>(ds);
  vtkPointSet::CopyStructure(ds);

  for (int i = 0; i < 3; i++)
    {
    this->Dimensions[i] = sg->Dimensions[i];
    }
  this->SetExtent(sg->GetExtent());
}

//----------------------------------------------------------------------------
void vtkRGrid::Initialize()
{
  this->Superclass::Initialize();

  if (this->Information)
    {
    this->SetDimensions(0, 0, 0);
    }
}

//----------------------------------------------------------------------------
int vtkRGrid::GetCellType(vtkIdType vtkNotUsed(cellId))
{
  return VTK_HEXAHEDRON;
}

//----------------------------------------------------------------------------
vtkIdType vtkRGrid::GetNumberOfCells()
{
  return this->Dimensions[0] * this->Dimensions[1] * this->Dimensions[2];
}

//----------------------------------------------------------------------------
void vtkRGrid::BuildLinks()
{
  // Remove the old links if they are already built
  if (this->Links)
    {
    this->Links->UnRegister(this);
    }

  this->Links = vtkCellLinks::New();
  this->Links->Allocate(this->GetNumberOfPoints());
  this->Links->Register(this);
  this->Links->BuildLinks(this, this->Cells);
  this->Links->Delete();
}

//----------------------------------------------------------------------------
vtkCell *vtkRGrid::GetCell(vtkIdType cellId)
{
  vtkCell* cell = this->Hexahedron;
  this->GetCell(cellId, static_cast<vtkCell*>(cell));
  return cell;
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCell(vtkIdType cellId, vtkGenericCell *cell)
{
  cell->SetCellTypeToHexahedron();
  this->GetCell(cellId, static_cast<vtkCell*>(cell));
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCell(vtkIdType cellId, vtkCell *cell)
{
  // Make sure data is defined
  if (!this->Points || !this->Cells)
    {
    vtkErrorMacro (<<"No data");
    }

  // Update dimensions
  this->GetDimensions();

  // Extract point coordinates and point ids. NOTE: the ordering of the
  // vtkHexahedron cells are tricky.
  vtkIdType *indices = this->GetCellPoints(cellId);
  for (int i = 0; i < 8; i++)
    {
    vtkIdType idx = indices[i];
    double x[3];
    this->Points->GetPoint(idx, x);
    cell->Points->SetPoint(i, x);
    cell->PointIds->SetId(i, idx);
    }
}

//----------------------------------------------------------------------------
// Fast implementation of GetCellBounds().  Bounds are calculated without
// constructing a cell.
void vtkRGrid::GetCellBounds(vtkIdType cellId, double bounds[6])
{
  // Make sure data is defined
  if (!this->Points)
    {
    vtkErrorMacro (<<"No data");
    return;
    }

  vtkMath::UninitializeBounds(bounds);

  // Update dimensions
  this->GetDimensions();

  vtkIdType *indices = this->GetCellPoints(cellId);
  double x[3];

  this->Points->GetPoint(indices[0], x);
  bounds[0] = bounds[1] = x[0];
  bounds[2] = bounds[3] = x[1];
  bounds[4] = bounds[5] = x[2];

  this->Points->GetPoint(indices[1], x);
  vtkAdjustBoundsMacro(bounds, x);

  this->Points->GetPoint(indices[2], x);
  vtkAdjustBoundsMacro(bounds, x);

  this->Points->GetPoint(indices[3], x);
  vtkAdjustBoundsMacro(bounds, x);

  this->Points->GetPoint(indices[4], x);
  vtkAdjustBoundsMacro(bounds, x);

  this->Points->GetPoint(indices[5], x);
  vtkAdjustBoundsMacro(bounds, x);

  this->Points->GetPoint(indices[6], x);
  vtkAdjustBoundsMacro(bounds, x);

  this->Points->GetPoint(indices[7], x);
  vtkAdjustBoundsMacro(bounds, x);
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCellDims(int cellDims[3])
{
  for (int i = 0; i < 3; ++i)
    {
    cellDims[i] = ((this->Dimensions[i]-1) < 1) ? 1 : this->Dimensions[i]-1;
    }
}

//----------------------------------------------------------------------------
// Set dimensions of structured grid dataset.
void vtkRGrid::SetDimensions(int i, int j, int k)
{
  this->SetExtent(0, i-1, 0, j-1, 0, k-1);
}

//----------------------------------------------------------------------------
// Set dimensions of structured grid dataset.
void vtkRGrid::SetDimensions(int dim[3])
{
  this->SetExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
}

//----------------------------------------------------------------------------
void vtkRGrid::GetPointCells(vtkIdType ptId, vtkIdList *cellIds)
{
  if (!this->Links)
    {
    this->BuildLinks();
    }
  cellIds->Reset();

  int numCells = this->Links->GetNcells(ptId);
  vtkIdType* cells = this->Links->GetCells(ptId);

  cellIds->SetNumberOfIds(numCells);
  for (int i = 0; i < numCells; i++)
    {
    cellIds->SetId(i, cells[i]);
    }
}

//----------------------------------------------------------------------------
vtkIdType* vtkRGrid::GetCellPoints(vtkIdType cellId)
{
  return this->Cells->GetData()->GetPointer(cellId * 9) + 1;
}

//----------------------------------------------------------------------------
// Get the points defining a cell. (See vtkDataSet for more info.)
void vtkRGrid::GetCellPoints(vtkIdType cellId, vtkIdList *ptIds)
{
  // Update dimensions
  this->GetDimensions();

  ptIds->Reset();

  ptIds->SetNumberOfIds(8);
  vtkIdType *indices = this->GetCellPoints(cellId);

  ptIds->SetId(0, indices[0]);
  ptIds->SetId(1, indices[1]);
  ptIds->SetId(2, indices[2]);
  ptIds->SetId(3, indices[3]);
  ptIds->SetId(4, indices[4]);
  ptIds->SetId(5, indices[5]);
  ptIds->SetId(6, indices[6]);
  ptIds->SetId(7, indices[7]);
}

//----------------------------------------------------------------------------
// Return a pointer to a list of point ids defining cell.
// More efficient than alternative method.
void vtkRGrid::GetCellPoints(vtkIdType cellId, vtkIdType& npts,
                             vtkIdType* &pts)
{
  pts = this->GetCellPoints(cellId);
  npts = 8;
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCellCoordinates(vtkIdType cellId, int &i, int &j, int &k)
{
  const int si = this->Dimensions[0];
  const int sj = this->Dimensions[1];
  i = j = k = 0;
  if (sj == 0)
    {
    i = cellId;
    }
  else if (si > 0 && sj > 0)
    {
    i = cellId / (si * sj);
    int jj = cellId % (si * sj);
    j = jj / si;
    k = jj % si;
    }
}

//----------------------------------------------------------------------------
void vtkRGrid::SetExtent(int extent[6])
{
  this->Modified();
  this->Dimensions[0] = extent[1] - extent[0] + 1;
  this->Dimensions[1] = extent[3] - extent[2] + 1;
  this->Dimensions[2] = extent[5] - extent[4] + 1;
  memcpy(this->Extent, extent, 6 * sizeof(int));
}

//----------------------------------------------------------------------------
void vtkRGrid::SetExtent(int xMin, int xMax,
                         int yMin, int yMax,
                         int zMin, int zMax)
{
  int extent[6];
  extent[0] = xMin; extent[1] = xMax;
  extent[2] = yMin; extent[3] = yMax;
  extent[4] = zMin; extent[5] = zMax;
  this->SetExtent(extent);
}

//----------------------------------------------------------------------------
int* vtkRGrid::GetDimensions()
{
  this->GetDimensions(this->Dimensions);
  return this->Dimensions;
}

//----------------------------------------------------------------------------
void vtkRGrid::GetDimensions(int dim[3])
{
  const int* extent = this->Extent;
  dim[0] = extent[1] - extent[0] + 1;
  dim[1] = extent[3] - extent[2] + 1;
  dim[2] = extent[5] - extent[4] + 1;
}

//----------------------------------------------------------------------------
// Determine neighbors as follows. Find the (shortest) list of cells that
// uses one of the points in ptIds. For each cell, in the list, see whether
// it contains the other points in the ptIds list. If so, it's a neighbor.
void vtkRGrid::GetCellNeighbors(vtkIdType cellId, vtkIdList *ptIds,
                                vtkIdList *cellIds)
{
  if (!this->Links)
    {
    this->BuildLinks();
    }

  cellIds->Reset();

  vtkIdType *minCells = 0;
  vtkIdType minPtId = 0;

  //Find the point used by the fewest number of cells

  vtkIdType numPts = ptIds->GetNumberOfIds();
  vtkIdType* pts = ptIds->GetPointer(0);
  vtkIdType minNumCells = VTK_INT_MAX;
  for (int i = 0; i < numPts; i++)
    {
    vtkIdType ptId = pts[i];
    vtkIdType numCells = this->Links->GetNcells(ptId);
    vtkIdType* cells = this->Links->GetCells(ptId);
    if ( numCells < minNumCells )
      {
      minNumCells = numCells;
      minCells = cells;
      minPtId = ptId;
      }
    }

  if (minNumCells == VTK_INT_MAX && numPts == 0)
    {
    vtkErrorMacro("input point ids empty.");
    minNumCells = 0;
    }

  //Now for each cell, see if it contains all the points
  //in the ptIds list.
  for (int i = 0; i < minNumCells; i++)
    {
    if ( minCells[i] != cellId ) //don't include current cell
      {
      vtkIdType* cellPts = this->GetCellPoints(minCells[i]);
      vtkIdType match = 1;
      for (int j = 0; j < numPts && match; j++) //for all pts in input cell
        {
        if (pts[j] != minPtId) //of course minPtId is contained by cell
          {
          match = 0;
          for (int k = 0; k < 8; k++) //for all points in candidate cell
            {
            if (pts[j] == cellPts[k])
              {
              match = 1; //a match was found
              break;
              }
            }//for all points in current cell
          }//if not guaranteed match
        }//for all points in input cell
      if (match)
        {
        cellIds->InsertNextId(minCells[i]);
        }
      }//if not the reference cell
    }//for all candidate cells attached to point
}

//----------------------------------------------------------------------------
unsigned long vtkRGrid::GetActualMemorySize()
{
  return this->vtkPointSet::GetActualMemorySize();
}

//----------------------------------------------------------------------------
void vtkRGrid::SetCells(vtkCellArray *cells)
{
  if (this->Cells == cells)
    {
    return;
    }
  if (this->Cells)
    {
    this->Cells->UnRegister(this);
    this->Cells = 0;
    }
  this->Cells = cells;
  if (this->Cells)
    {
    this->Cells->Register(this);
    }
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkRGrid::ShallowCopy(vtkDataObject *dataObject)
{
  vtkRGrid *grid = vtkRGrid::SafeDownCast(dataObject);

  if (grid)
    {
    this->InternalRGridCopy(grid);
    this->CellVisibility->ShallowCopy(grid->CellVisibility);
    }

  // Do superclass
  this->vtkPointSet::ShallowCopy(dataObject);
}

//----------------------------------------------------------------------------
void vtkRGrid::DeepCopy(vtkDataObject *dataObject)
{
  vtkRGrid* grid = vtkRGrid::SafeDownCast(dataObject);

  if (grid)
    {
    this->InternalRGridCopy(grid);
    this->CellVisibility->DeepCopy(grid->CellVisibility);
    }

  // Do superclass
  this->vtkPointSet::DeepCopy(dataObject);
}

//----------------------------------------------------------------------------
// This copies all the local variables (but not objects).
void vtkRGrid::InternalRGridCopy(vtkRGrid *src)
{
  // Update dimensions
  this->GetDimensions();

  memcpy(this->Dimensions, src->Dimensions, 3 * sizeof(int));
  memcpy(this->Extent, src->GetExtent(), 6 * sizeof(int));
}

//----------------------------------------------------------------------------
// Override this method because of blanking
void vtkRGrid::ComputeScalarRange()
{
  if (this->GetMTime() > this->ScalarRangeComputeTime)
    {
    vtkDataArray *ptScalars = this->PointData->GetScalars();
    vtkDataArray *cellScalars = this->CellData->GetScalars();

    double ptRange[2];
    ptRange[0] = VTK_DOUBLE_MAX;
    ptRange[1] = VTK_DOUBLE_MIN;
    if (ptScalars)
      {
      int num = this->GetNumberOfPoints();
      for (int id = 0; id < num; id++)
        {
        double s = ptScalars->GetComponent(id, 0);
        if (s < ptRange[0])
          {
          ptRange[0] = s;
          }
        if (s > ptRange[1])
          {
          ptRange[1] = s;
          }
        }
      }

    double cellRange[2];
    cellRange[0] = ptRange[0];
    cellRange[1] = ptRange[1];
    if (cellScalars)
      {
      int num = this->GetNumberOfCells();
      for (int id = 0; id < num; id++)
        {
        double s = cellScalars->GetComponent(id, 0);
        if (s < cellRange[0])
          {
          cellRange[0] = s;
          }
        if (s > cellRange[1])
          {
          cellRange[1] = s;
          }
        }
      }

    this->ScalarRange[0] = (cellRange[0] >= VTK_DOUBLE_MAX ? 0.0 : cellRange[0]);
    this->ScalarRange[1] = (cellRange[1] <= VTK_DOUBLE_MIN ? 1.0 : cellRange[1]);

    this->ScalarRangeComputeTime.Modified();
    }
}

//----------------------------------------------------------------------------
// Turn off a particular data cell.
void vtkRGrid::BlankCell(vtkIdType cellId)
{
  int celldims[3];
  this->GetCellDims(celldims);
  this->CellVisibility->Initialize(celldims);
  this->CellVisibility->Blank(cellId);
}

//----------------------------------------------------------------------------
// Turn on a particular data cell.
void vtkRGrid::UnBlankCell(vtkIdType cellId)
{
  int celldims[3];
  this->GetCellDims(celldims);
  this->CellVisibility->Initialize(celldims);
  this->CellVisibility->UnBlank(cellId);
}

//----------------------------------------------------------------------------
void vtkRGrid::SetCellVisibilityArray(vtkUnsignedCharArray *cellVis)
{
  this->CellVisibility->SetVisibilityById(cellVis);
}

//----------------------------------------------------------------------------
vtkUnsignedCharArray* vtkRGrid::GetCellVisibilityArray()
{
  int celldims[3];
  this->GetCellDims(celldims);
  this->CellVisibility->Initialize(celldims);
  this->CellVisibility->Allocate();
  return this->CellVisibility->GetVisibilityById();
}

//----------------------------------------------------------------------------
// Return non-zero if the specified cell is visible (i.e., not blanked)
unsigned char vtkRGrid::IsCellVisible(vtkIdType cellId)
{
  return this->CellVisibility->IsVisible(cellId);
}

//----------------------------------------------------------------------------
unsigned char vtkRGrid::GetCellBlanking()
{
  return this->CellVisibility->IsConstrained();
}

//----------------------------------------------------------------------------
void vtkRGrid::Crop(const int* updateExtent)
{
  int uExt[6];
  const int* extent = this->Extent;

  // If the update extent is larger than the extent,
  // we cannot do anything about it here.
  for (int i = 0; i < 3; ++i)
    {
    uExt[i*2] = updateExtent[i*2];
    if (uExt[i*2] < extent[i*2])
      {
      uExt[i*2] = extent[i*2];
      }
    uExt[i*2+1] = updateExtent[i*2+1];
    if (uExt[i*2+1] > extent[i*2+1])
      {
      uExt[i*2+1] = extent[i*2+1];
      }
    }

  // If extents already match, then we need to do nothing.
  if (extent[0] == uExt[0] && extent[1] == uExt[1]
      && extent[2] == uExt[2] && extent[3] == uExt[3]
      && extent[4] == uExt[4] && extent[5] == uExt[5])
    {
    return;
    }
  else
    {
    // Get the points.  Protect against empty data objects.
    vtkPoints *inPts = this->GetPoints();
    if (!inPts)
      {
      return;
      }

    vtkDebugMacro(<< "Cropping RGrid");

    vtkCellData*  inCD  = this->GetCellData();
    vtkNew<vtkCellData> outCD;

    // Allocate necessary objects
    int outSize = (uExt[1] - uExt[0] + 1) *
      (uExt[3] - uExt[2] + 1) *
      (uExt[5] - uExt[4] + 1);
    outCD->CopyAllocate(inCD, outSize, outSize);

    vtkNew<vtkCellArray> cells;
    cells->Allocate(outSize * 9);

    // Traverse input data and copy cell attributes to output
    int inInc1 = (extent[1] - extent[0]);
    int inInc2 = inInc1 * (extent[3] - extent[2]);
    for (int k = uExt[4]; k < uExt[5]; ++k)
      {
      int kOffset = (k - extent[4]) * inInc2;
      for (int j = uExt[2]; j < uExt[3]; ++j)
        {
        int jOffset = (j - extent[2]) * inInc1;
        for (int i = uExt[0]; i < uExt[1]; ++i)
          {
          int idx = (i - extent[0]) + jOffset + kOffset;
          vtkNew<vtkIdList> ptIds;
          this->GetCellPoints(idx, ptIds.GetPointer());
          vtkIdType nCellId = cells->InsertNextCell(ptIds.GetPointer());
          outCD->CopyData(inCD, idx, nCellId);
          }
        }
      }

    this->SetExtent(uExt);
    inCD->ShallowCopy(outCD.GetPointer());
    }
}

//----------------------------------------------------------------------------
void vtkRGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  int dim[3];
  this->GetDimensions(dim);
  os << indent << "Dimensions: (" << dim[0] << ", "
                                  << dim[1] << ", "
                                  << dim[2] << ")\n";

  const int* extent = this->Extent;
  os << indent << "Extent: "
    << extent[0] << ", " << extent[1] << ", "
    << extent[2] << ", " << extent[3] << ", "
    << extent[4] << ", " << extent[5] << endl;

  os << ")\n";
}

//----------------------------------------------------------------------------
vtkRGrid* vtkRGrid::GetData(vtkInformation* info)
{
  return info? vtkRGrid::SafeDownCast(info->Get(DATA_OBJECT())) : 0;
}

//----------------------------------------------------------------------------
vtkRGrid* vtkRGrid::GetData(vtkInformationVector* v, int i)
{
  return vtkRGrid::GetData(v->GetInformationObject(i));
}
