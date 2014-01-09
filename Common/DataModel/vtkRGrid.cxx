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

#include "vtkCellData.h"
#include "vtkEmptyCell.h"
#include "vtkGenericCell.h"
#include "vtkHexahedron.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkQuad.h"
#include "vtkVertex.h"

vtkStandardNewMacro(vtkRGrid);

#define vtkAdjustBoundsMacro( A, B ) \
  A[0] = (B[0] < A[0] ? B[0] : A[0]);   A[1] = (B[0] > A[1] ? B[0] : A[1]); \
  A[2] = (B[1] < A[2] ? B[1] : A[2]);   A[3] = (B[1] > A[3] ? B[1] : A[3]); \
  A[4] = (B[2] < A[4] ? B[2] : A[4]);   A[5] = (B[2] > A[5] ? B[2] : A[5])

vtkRGrid::vtkRGrid()
{
  this->Hexahedron = vtkHexahedron::New();
  this->EmptyCell = vtkEmptyCell::New();

  this->Dimensions[0] = 0;
  this->Dimensions[1] = 0;
  this->Dimensions[2] = 0;

  int extent[6] = {0, -1, 0, -1, 0, -1};
  memcpy(this->Extent, extent, 6*sizeof(int));

  this->Information->Set(vtkDataObject::DATA_EXTENT_TYPE(), VTK_3D_EXTENT);
  this->Information->Set(vtkDataObject::DATA_EXTENT(), this->Extent, 6);
}

//----------------------------------------------------------------------------
vtkRGrid::~vtkRGrid()
{
  this->Hexahedron->Delete();
  this->EmptyCell->Delete();
}

//----------------------------------------------------------------------------
// Copy the geometric and topological structure of an input structured grid.
void vtkRGrid::CopyStructure(vtkDataSet *ds)
{
  vtkRGrid *sg=static_cast<vtkRGrid *>(ds);
  vtkPointSet::CopyStructure(ds);
  int i;

  for (i=0; i<3; i++)
    {
    this->Dimensions[i] = sg->Dimensions[i];
    }
  this->SetExtent(sg->GetExtent());
}


//----------------------------------------------------------------------------
void vtkRGrid::Initialize()
{
  this->Superclass::Initialize();

  if(this->Information)
    {
    this->SetDimensions(0,0,0);
    }
}

//----------------------------------------------------------------------------
int vtkRGrid::GetCellType(vtkIdType vtkNotUsed(cellId))
{
  return VTK_HEXAHEDRON;
}

//----------------------------------------------------------------------------
vtkCell *vtkRGrid::GetCell(vtkIdType cellId)
{
  vtkCell *cell = NULL;
  vtkIdType idx;
  int i, j, k;
  int d01, offset1, offset2;

  // Make sure data is defined
  if ( ! this->Points )
    {
    vtkErrorMacro (<<"No data");
    return NULL;
    }

  // Update dimensions
  this->GetDimensions();

  cell = this->Hexahedron;
  d01 = this->Dimensions[0]*this->Dimensions[1];
  i = cellId % (this->Dimensions[0] - 1);
  j = (cellId / (this->Dimensions[0] - 1)) % (this->Dimensions[1] - 1);
  k = cellId / ((this->Dimensions[0] - 1) * (this->Dimensions[1] - 1));
  idx = i+ j*this->Dimensions[0] + k*d01;
  offset1 = 1;
  offset2 = this->Dimensions[0];

  cell->PointIds->SetId(0,idx);
  cell->PointIds->SetId(1,idx+offset1);
  cell->PointIds->SetId(2,idx+offset1+offset2);
  cell->PointIds->SetId(3,idx+offset2);
  idx += d01;
  cell->PointIds->SetId(4,idx);
  cell->PointIds->SetId(5,idx+offset1);
  cell->PointIds->SetId(6,idx+offset1+offset2);
  cell->PointIds->SetId(7,idx+offset2);

  // Extract point coordinates and point ids. NOTE: the ordering of the vtkQuad
  // and vtkHexahedron cells are tricky.
  int NumberOfIds = cell->PointIds->GetNumberOfIds();
  for (i=0; i<NumberOfIds; i++)
    {
    idx = cell->PointIds->GetId(i);
    cell->Points->SetPoint(i,this->Points->GetPoint(idx));
    }

  return cell;
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCell(vtkIdType cellId, vtkGenericCell *cell)
{
  vtkIdType   idx;
  int   i, j, k;
  int   d01, offset1, offset2;
  double x[3];

  // Make sure data is defined
  if ( ! this->Points )
    {
    vtkErrorMacro (<<"No data");
    }

  // Update dimensions
  this->GetDimensions();

  cell->SetCellTypeToHexahedron();
  d01 = this->Dimensions[0]*this->Dimensions[1];
  i = cellId % (this->Dimensions[0] - 1);
  j = (cellId / (this->Dimensions[0] - 1)) % (this->Dimensions[1] - 1);
  k = cellId / ((this->Dimensions[0] - 1) * (this->Dimensions[1] - 1));
  idx = i+ j*this->Dimensions[0] + k*d01;
  offset1 = 1;
  offset2 = this->Dimensions[0];

  cell->PointIds->SetId(0,idx);
  cell->PointIds->SetId(1,idx+offset1);
  cell->PointIds->SetId(2,idx+offset1+offset2);
  cell->PointIds->SetId(3,idx+offset2);
  idx += d01;
  cell->PointIds->SetId(4,idx);
  cell->PointIds->SetId(5,idx+offset1);
  cell->PointIds->SetId(6,idx+offset1+offset2);
  cell->PointIds->SetId(7,idx+offset2);

  // Extract point coordinates and point ids. NOTE: the ordering of the vtkQuad
  // and vtkHexahedron cells are tricky.
  int NumberOfIds = cell->PointIds->GetNumberOfIds();
  for (i=0; i<NumberOfIds; i++)
    {
    idx = cell->PointIds->GetId(i);
    this->Points->GetPoint(idx, x);
    cell->Points->SetPoint(i, x);
    }
}


//----------------------------------------------------------------------------
// Fast implementation of GetCellBounds().  Bounds are calculated without
// constructing a cell.
void vtkRGrid::GetCellBounds(vtkIdType cellId, double bounds[6])
{
  vtkIdType idx = 0;
  int i, j, k;
  vtkIdType d01;
  int offset1 = 0;
  int offset2 = 0;
  double x[3];

  // Make sure data is defined
  if ( ! this->Points )
    {
    vtkErrorMacro (<<"No data");
    return;
    }

  vtkMath::UninitializeBounds(bounds);

  // Update dimensions
  this->GetDimensions();

  d01 = this->Dimensions[0]*this->Dimensions[1];
  i = cellId % (this->Dimensions[0] - 1);
  j = (cellId / (this->Dimensions[0] - 1)) % (this->Dimensions[1] - 1);
  k = cellId / ((this->Dimensions[0] - 1) * (this->Dimensions[1] - 1));
  idx = i+ j*this->Dimensions[0] + k*d01;
  offset1 = 1;
  offset2 = this->Dimensions[0];

  this->Points->GetPoint(idx, x);
  bounds[0] = bounds[1] = x[0];
  bounds[2] = bounds[3] = x[1];
  bounds[4] = bounds[5] = x[2];

  this->Points->GetPoint( idx+offset1, x);
  vtkAdjustBoundsMacro( bounds, x );

  this->Points->GetPoint( idx+offset1+offset2, x);
  vtkAdjustBoundsMacro( bounds, x );

  this->Points->GetPoint( idx+offset2, x);
  vtkAdjustBoundsMacro( bounds, x );

  idx += d01;

  this->Points->GetPoint(idx, x);
  vtkAdjustBoundsMacro( bounds, x );

  this->Points->GetPoint( idx+offset1, x);
  vtkAdjustBoundsMacro( bounds, x );

  this->Points->GetPoint( idx+offset1+offset2, x);
  vtkAdjustBoundsMacro( bounds, x );

  this->Points->GetPoint( idx+offset2, x);
  vtkAdjustBoundsMacro( bounds, x );
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCellDims( int cellDims[3] )
{
  for( int i=0; i < 3; ++i )
    {
    cellDims[i] = ( (this->Dimensions[i]-1) < 1)? 1 : this->Dimensions[i]-1;
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
}

//----------------------------------------------------------------------------
// Get the points defining a cell. (See vtkDataSet for more info.)
void vtkRGrid::GetCellPoints(vtkIdType cellId, vtkIdList *ptIds)
{
  // Update dimensions
  this->GetDimensions();

  int iMin, iMax, jMin, jMax, kMin, kMax;
  vtkIdType d01 = this->Dimensions[0]*this->Dimensions[1];

  ptIds->Reset();
  iMin = iMax = jMin = jMax = kMin = kMax = 0;

  iMin = cellId % (this->Dimensions[0] - 1);
  iMax = iMin + 1;
  jMin = (cellId / (this->Dimensions[0] - 1)) % (this->Dimensions[1] - 1);
  jMax = jMin + 1;
  kMin = cellId / ((this->Dimensions[0] - 1) * (this->Dimensions[1] - 1));
  kMax = kMin + 1;
  ptIds->SetNumberOfIds(8);
  ptIds->SetId(0, iMin + jMin*this->Dimensions[0] + kMin*d01);
  ptIds->SetId(1, iMax + jMin*this->Dimensions[0] + kMin*d01);
  ptIds->SetId(2, iMax + jMax*this->Dimensions[0] + kMin*d01);
  ptIds->SetId(3, iMin + jMax*this->Dimensions[0] + kMin*d01);
  ptIds->SetId(4, iMin + jMin*this->Dimensions[0] + kMax*d01);
  ptIds->SetId(5, iMax + jMin*this->Dimensions[0] + kMax*d01);
  ptIds->SetId(6, iMax + jMax*this->Dimensions[0] + kMax*d01);
  ptIds->SetId(7, iMin + jMax*this->Dimensions[0] + kMax*d01);
}

//----------------------------------------------------------------------------
void vtkRGrid::SetExtent(int extent[6])
{
  this->Modified();
  this->Dimensions[0] = extent[1] - extent[0] + 1;
  this->Dimensions[1] = extent[3] - extent[2] + 1;
  this->Dimensions[2] = extent[5] - extent[4] + 1;
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

int *vtkRGrid::GetDimensions ()
{
  this->GetDimensions(this->Dimensions);
  return this->Dimensions;
}

void vtkRGrid::GetDimensions (int dim[3])
{
  const int* extent = this->Extent;
  dim[0] = extent[1] - extent[0] + 1;
  dim[1] = extent[3] - extent[2] + 1;
  dim[2] = extent[5] - extent[4] + 1;
}

//----------------------------------------------------------------------------
void vtkRGrid::GetCellNeighbors(vtkIdType cellId, vtkIdList *ptIds,
                                         vtkIdList *cellIds)
{
  int numPtIds=ptIds->GetNumberOfIds();

  // Use special methods for speed
  switch (numPtIds)
    {
    case 0:
      cellIds->Reset();
      return;

    case 1: case 2: case 4: //vertex, edge, face neighbors
      //vtkRData::GetCellNeighbors(cellId, ptIds,
                                          //cellIds, this->GetDimensions());
      break;

    default:
      this->vtkDataSet::GetCellNeighbors(cellId, ptIds, cellIds);
    }
}

//----------------------------------------------------------------------------
unsigned long vtkRGrid::GetActualMemorySize()
{
  return this->vtkPointSet::GetActualMemorySize();
}

//----------------------------------------------------------------------------
void vtkRGrid::ShallowCopy(vtkDataObject *dataObject)
{
  vtkRGrid *grid = vtkRGrid::SafeDownCast(dataObject);

  if ( grid != NULL )
    {
    this->InternalRGridCopy(grid);
    }


  // Do superclass
  this->vtkPointSet::ShallowCopy(dataObject);
}

//----------------------------------------------------------------------------
void vtkRGrid::DeepCopy(vtkDataObject *dataObject)
{
  vtkRGrid *grid = vtkRGrid::SafeDownCast(dataObject);

  if ( grid != NULL )
    {
    this->InternalRGridCopy(grid);
    }

  // Do superclass
  this->vtkPointSet::DeepCopy(dataObject);
}

//----------------------------------------------------------------------------
// This copies all the local variables (but not objects).
void vtkRGrid::InternalRGridCopy(vtkRGrid *src)
{
  int idx;

  // Update dimensions
  this->GetDimensions();

  for (idx = 0; idx < 3; ++idx)
    {
    this->Dimensions[idx] = src->Dimensions[idx];
    }
  memcpy(this->Extent, src->GetExtent(), 6*sizeof(int));
}

//----------------------------------------------------------------------------
// Override this method because of blanking
void vtkRGrid::ComputeScalarRange()
{
  if ( this->GetMTime() > this->ScalarRangeComputeTime )
    {
    vtkDataArray *ptScalars = this->PointData->GetScalars();
    vtkDataArray *cellScalars = this->CellData->GetScalars();
    double ptRange[2];
    double cellRange[2];
    double s;
    int id, num;

    ptRange[0] =  VTK_DOUBLE_MAX;
    ptRange[1] =  VTK_DOUBLE_MIN;
    if ( ptScalars )
      {
      num = this->GetNumberOfPoints();
      for (id=0; id < num; id++)
        {
        s = ptScalars->GetComponent(id,0);
        if ( s < ptRange[0] )
          {
          ptRange[0] = s;
          }
        if ( s > ptRange[1] )
          {
          ptRange[1] = s;
          }
        }
      }

    cellRange[0] =  ptRange[0];
    cellRange[1] =  ptRange[1];
    if ( cellScalars )
      {
      num = this->GetNumberOfCells();
      for (id=0; id < num; id++)
        {
        s = cellScalars->GetComponent(id,0);
        if ( s < cellRange[0] )
          {
          cellRange[0] = s;
          }
        if ( s > cellRange[1] )
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
void vtkRGrid::Crop(const int* updateExtent)
{
  int i, j, k;
  int uExt[6];
  const int* extent = this->Extent;

  // If the update extent is larger than the extent,
  // we cannot do anything about it here.
  for (i = 0; i < 3; ++i)
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
    vtkRGrid *newGrid;
    vtkPointData *inPD, *outPD;
    vtkCellData *inCD, *outCD;
    int outSize, jOffset, kOffset;
    vtkIdType idx, newId;
    vtkPoints *newPts, *inPts;
    int inInc1, inInc2;

    // Get the points.  Protect against empty data objects.
    inPts = this->GetPoints();
    if (inPts == NULL)
      {
      return;
      }

    vtkDebugMacro(<< "Cropping Grid");

    newGrid = vtkRGrid::New();
    inPD  = this->GetPointData();
    inCD  = this->GetCellData();
    outPD = newGrid->GetPointData();
    outCD = newGrid->GetCellData();

    // Allocate necessary objects
    //
    newGrid->SetExtent(uExt);
    outSize = (uExt[1]-uExt[0]+1)*(uExt[3]-uExt[2]+1)*(uExt[5]-uExt[4]+1);
    newPts = inPts->NewInstance();
    newPts->SetDataType(inPts->GetDataType());
    newPts->SetNumberOfPoints(outSize);
    outPD->CopyAllocate(inPD,outSize,outSize);
    outCD->CopyAllocate(inCD,outSize,outSize);

    // Traverse this data and copy point attributes to output
    newId = 0;
    inInc1 = (extent[1]-extent[0]+1);
    inInc2 = inInc1*(extent[3]-extent[2]+1);
    for ( k=uExt[4]; k <= uExt[5]; ++k)
      {
      kOffset = (k - extent[4]) * inInc2;
      for ( j=uExt[2]; j <= uExt[3]; ++j)
        {
        jOffset = (j - extent[2]) * inInc1;
        for ( i=uExt[0]; i <= uExt[1]; ++i)
          {
          idx = (i - extent[0]) + jOffset + kOffset;
          newPts->SetPoint(newId,inPts->GetPoint(idx));
          outPD->CopyData(inPD, idx, newId++);
          }
        }
      }

    // Traverse input data and copy cell attributes to output
    newId = 0;
    inInc1 = (extent[1] - extent[0]);
    inInc2 = inInc1*(extent[3] - extent[2]);
    for ( k=uExt[4]; k < uExt[5]; ++k )
      {
      kOffset = (k - extent[4]) * inInc2;
      for ( j=uExt[2]; j < uExt[3]; ++j )
        {
        jOffset = (j - extent[2]) * inInc1;
        for ( i=uExt[0]; i < uExt[1]; ++i )
          {
          idx = (i - extent[0]) + jOffset + kOffset;
          outCD->CopyData(inCD, idx, newId++);
          }
        }
      }

    this->SetExtent(uExt);
    this->SetPoints(newPts);
    newPts->Delete();
    inPD->ShallowCopy(outPD);
    inCD->ShallowCopy(outCD);
    newGrid->Delete();
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
  os << indent << "Extent: " << extent[0] << ", "
     << extent[1] << ", " << extent[2] << ", "
     << extent[3] << ", " << extent[4] << ", "
     << extent[5] << endl;

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

//----------------------------------------------------------------------------
void vtkRGrid::GetPoint(
    int i, int j, int k, double p[3], bool adjustForExtent)
{
  int extent[6];
  this->GetExtent(extent);

  if(i < extent[0] || i > extent[1] ||
     j < extent[2] || j > extent[3] ||
     k < extent[4] || k > extent[5])
    {
    vtkErrorMacro("ERROR: IJK coordinates are outside of grid extent!");
    return; // out of bounds!
    }

  int pos[3];
  pos[0] = i;
  pos[1] = j;
  pos[2] = k;

  vtkIdType id = 0;

  if(adjustForExtent)
    {
    //id = vtkRData::ComputePointIdForExtent(extent, pos);
    }
  else
    {
    int dim[3];
    this->GetDimensions(dim);
    //id = vtkRData::ComputePointId(dim, pos);
    }

  this->GetPoint(id, p);
}
