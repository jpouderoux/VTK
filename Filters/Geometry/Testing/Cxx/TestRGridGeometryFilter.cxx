/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestUnstructuredGridGeometryFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkLookupTable.h"
#include "vtkPolyData.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFloatArray.h"
#include "vtkGeometryFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkCamera.h"
#include "vtkPointData.h"
#include "vtkProperty.h"
#include "vtkRGrid.h"
#include "vtkRGridSlice.h"
#include "vtkPoints.h"
#include "vtkShrinkFilter.h"
#include "vtkCell.h"
#include "vtkNew.h"
#include "vtkOutlineCornerFilter.h"
#include "vtkInteractorStyleTrackballCamera.h"
//#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkObjectFactory.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkPointLocator.h"
#include "vtkMergePoints.h"
#include "vtkCleanPolyData.h"
#include "vtkLODProp3D.h"
#include "vtkDecimatePro.h"
#include "vtkOutlineFilter.h"
#include "vtkLODActor.h"
#include "vtkTriangleFilter.h"
#include "vtkStripper.h"
#include "vtkQuadricClustering.h"
#include "vtkQuadricLODActor.h"

vtkNew<vtkRenderer> renderer;
vtkNew<vtkRenderWindow> renWin;
vtkNew<vtkRGridSlice> slice;

//----------------------------------------------------------------------------
void CleanUnstructuredGrid(vtkUnstructuredGrid *input)
{
  if (input->GetNumberOfCells() == 0)
    {
    return;
    }

  // First, create a new points array that eliminate duplicate points.
  // Also create a mapping from the old point id to the new.
  vtkNew<vtkPoints> newPts;
  vtkIdType num = input->GetNumberOfPoints();
  vtkIdType* ptMap = new vtkIdType[num];

  vtkNew<vtkPointLocator> locator;
  locator->InitPointInsertion(newPts.GetPointer(), input->GetBounds(), num);

  for (vtkIdType id = 0; id < num; ++id)
    {
    double pt[3];
    input->GetPoint(id, pt);
    vtkIdType newId;
    if (locator->InsertUniquePoint(pt, newId))
      {
      //output->GetPointData()->CopyData(input->GetPointData(), id, newId);
      }
    ptMap[id] = newId;
    }
  input->SetPoints(newPts.GetPointer());

  // Now copy the cells.
  vtkNew<vtkIdList> cellPoints;
  num = input->GetNumberOfCells();
  vtkNew<vtkCellArray> newCells;
  newCells->Allocate(num * 9);
  for (vtkIdType id = 0; id < num; ++id)
    {
    input->GetCellPoints(id, cellPoints.GetPointer());
    for (int i = 0; i < cellPoints->GetNumberOfIds(); i++)
      {
      int cellPtId = cellPoints->GetId(i);
      cellPoints->SetId(i, ptMap[cellPtId]);
      }
    newCells->InsertNextCell(cellPoints.GetPointer());
    }

  delete [] ptMap;

  newCells->Squeeze();
  input->SetCells(VTK_HEXAHEDRON, newCells.GetPointer());
}

//-----------------------------------------------------------------------------
class MytInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static MytInteractorStyle* New();
  vtkTypeMacro(MytInteractorStyle, vtkInteractorStyleTrackballCamera);

  MytInteractorStyle()
  {
  }

  virtual void OnKeyPress()
  {
    vtkRenderWindowInteractor* rwi = this->GetInteractor();
    std::string key = rwi->GetKeySym();
    char keyCode = rwi->GetKeyCode();
    std::string keyname = rwi->GetKeySym();
    int extents[6];
    slice->GetSliceExtents(extents);
    if (keyname == "Left")
      {
      extents[0] -= 1;
      extents[1] = extents[0] + 1;
      extents[2] = 0;
      extents[3] = VTK_INT_MAX;
      }
    else if (keyname == "Right")
      {
      extents[0] += 1;
      extents[1] = extents[0] + 1;
      extents[2] = 0;
      extents[3] = VTK_INT_MAX;
      }
    if (keyname == "Down")
      {
      extents[0] = 0;
      extents[1] = VTK_INT_MAX;
      extents[2] -= 1;
      extents[3] = extents[2] + 1;
      }
    else if (keyname == "Up")
      {
      extents[0] = 0;
      extents[1] = VTK_INT_MAX;
      extents[2] += 1;
      extents[3] = extents[2] + 1;
      }
    else if (keyname == "Prior")
      {
      extents[0] = 0;
      extents[1] = VTK_INT_MAX;
      extents[4] += 1;
      extents[5] = extents[4] + 1;
      }
    else if (keyname == "Next")
      {
      extents[0] = 0;
      extents[1] = VTK_INT_MAX;
      extents[4] -= 1;
      extents[5] = extents[4] + 1;
      }
    /*else switch (keyCode)
      {
      default:break;
      }*/

    slice->SetSliceExtents(extents);
    slice->Update();
    renWin->Render();
  }
};

void CreateGrid(vtkRGrid* grid)
{
  vtkNew<vtkPoints> pts;
  const int si = 60;
  const int sj = 50;
  const int sk = 20;
  const int si1 = si+3;
  const int sj1 = sj+2;
  const int sk1 = sk+1;

  pts->SetNumberOfPoints(si1 * sj1 * sk1);
  for (int k = 0, ijk = 0; k < sk1; k++)
    {
    for (int j = 0; j < sj1; j++)
      {
      for (int i = 0; i < si1; i++, ijk++)
        {
        double x = i;
        double y = j;
        double z = k;
        if (i > 30)
          {
          x -= 1;
          z = k + (j*8. / 52.);//1.0;
          if (i > 50)
            {
            x -= 1;
            z = k + 2.4;
            }
          }
        if (j > 20)
          {
          y -= 1;
          z += 2.2 * (j*8. / 52.);
          }
        pts->SetPoint(ijk, x, y, z);
        }
      }
    }

  vtkNew<vtkDoubleArray> scalars;
  scalars->SetName("scalars");
  scalars->SetNumberOfValues(si * sj * sk);

  vtkNew<vtkUnsignedCharArray> blanking;
  blanking->SetName("blanking");
  blanking->SetNumberOfValues(si * sj * sk);

  vtkNew<vtkCellArray> cells;
  cells->Allocate(9 * si * sj * sk);

  for (int k = 0, ijk = 0; k < sk; k++)
    {
    for (int j = 0; j < sj; j++)
      {
      for (int i = 0; i < si; i++, ijk++)
        {
        int x = i;
        int y = j;
        int z = k;
        if (i >= 30)
          {
          x++;
          if (i >= 49)
            {
            x++;
            }
          }
        if (j >= 20)
          {
          y++;
          }
#define CTOI(_i,_j,_k) _i + (_j) * si1 + (_k) * sj1 * si1

        cells->InsertNextCell(8);
        for (int kk = 0; kk < 2; kk++)
          {
          cells->InsertCellPoint(CTOI(x,   y,   z+kk));
          cells->InsertCellPoint(CTOI(x+1, y,   z+kk));
          cells->InsertCellPoint(CTOI(x+1, y+1, z+kk));
          cells->InsertCellPoint(CTOI(x,   y+1, z+kk));
          }

        scalars->SetValue(ijk, k);
        bool visible = true;
        if (i > 4 && i < 18 && j < 10 && k < 10)
          visible = false;

        blanking->SetValue(ijk, visible ? 1 : 0);
        }
      }
    }

  grid->SetDimensions(si, sj, sk);
  grid->SetPoints(pts.GetPointer());
  grid->SetCells(cells.GetPointer());
  grid->GetCellData()->AddArray(blanking.GetPointer());
  grid->GetCellData()->AddArray(scalars.GetPointer());
  grid->GetCellData()->SetActiveScalars("scalars");
  grid->SetCellVisibilityArray(blanking.GetPointer());
}

void LoadGrid(vtkRGrid* grid2, const std::string& fname, double dataRange[2])
{
  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(fname.c_str());
  reader->Update();
  vtkUnstructuredGrid* ug = reader->GetOutput();

  CleanUnstructuredGrid(ug);

  grid2->SetPoints(ug->GetPoints());
  vtkIntArray* ijk =
    vtkIntArray::SafeDownCast(ug->GetCellData()->GetArray("IJK"));
  vtkCellArray* ugcells = ug->GetCells();
  ugcells->InitTraversal();

  vtkIdType npts, *apts, cid = 0;
  int dims[3] = { 0, 0, 0 };
  while(ugcells->GetNextCell(npts, apts))
    {
    double* ijk_ = ijk->GetTuple3(cid);
    if (dims[0] < ijk_[0]) dims[0] = ijk_[0];
    if (dims[1] < ijk_[1]) dims[1] = ijk_[1];
    if (dims[2] < ijk_[2]) dims[2] = ijk_[2];

    cid++;
    }
  dims[0]++; dims[1]++; dims[2]++;

  vtkNew<vtkUnsignedCharArray> blanking2;
  blanking2->SetName("blanking");
  vtkIdType len = ijk->GetNumberOfTuples();
  blanking2->SetNumberOfValues(len); //dims[0] * dims[1] * dims[2]);
  vtkFloatArray* por = vtkFloatArray::SafeDownCast(ug->GetCellData()->GetArray("PORV"));
  dataRange[0] = VTK_DOUBLE_MAX;
  dataRange[1] = VTK_DOUBLE_MIN;
  for (int i = 0; i < len; i++)
    {
    double v = por->GetValue(i);
    if (v != 0. && v < dataRange[0]) dataRange[0] = v;
    if (v != 0. && v > dataRange[1]) dataRange[1] = v;
    blanking2->SetValue(i, (v != 0.) ? 1 : 0);
    }

  grid2->SetDimensions(dims[0], dims[1], dims[2]);
  grid2->SetCells(ugcells);
  grid2->GetCellData()->ShallowCopy(ug->GetCellData());
  grid2->GetPointData()->ShallowCopy(ug->GetPointData());
  grid2->GetCellData()->AddArray(blanking2.GetPointer());
  grid2->SetCellVisibilityArray(blanking2.GetPointer());
  grid2->GetCellData()->SetActiveScalars("PORV");
}

vtkStandardNewMacro(MytInteractorStyle);

int TestRGridGeometryFilter(int argc, char* argv[])
{
  double dataRange[2];
  vtkNew<vtkRGrid> grid;
  //CreateGrid(grid.GetPointer());
  //LoadGrid(grid.GetPointer(), "c:/BRILLIG.bin.vtu", dataRange);
  LoadGrid(grid.GetPointer(), "c:/BEST_700M_60AC_FAULT.bin.vtu", dataRange);

  cout << "Grid has " << grid->GetNumberOfCells() << " cells and "
    << grid->GetNumberOfPoints() << " pts "<< endl;

  vtkPolyData* pd = 0;
    {
    vtkNew<vtkDataSetSurfaceFilter> surface;
    surface->SetInputData(grid.GetPointer());
    surface->Update();
    pd = surface->GetOutput();
    cout << "Geometry has " << pd->GetNumberOfCells() << " cells and "
      << pd->GetNumberOfPoints() << " points "<< endl;
    vtkNew<vtkCleanPolyData> clean;
    clean->SetInputData(pd);
    clean->Update();
    pd = clean->GetOutput();
    pd->Register(0);
    }

  cout << "Geometry has " << pd->GetNumberOfCells() << " cells and "
    << pd->GetNumberOfPoints() << " points "<< endl;
  pd->PrintSelf(cout, vtkIndent());

  slice->SetInputData(grid.GetPointer());
  vtkNew<vtkDataSetSurfaceFilter> slsurface;
  slsurface->SetInputConnection(slice->GetOutputPort(0));

  // Rendering pipeline
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(pd);
  mapper->SetScalarRange(dataRange);

  vtkNew<vtkQuadricLODActor> lodActor;
  lodActor->SetMapper(mapper.GetPointer());
  lodActor->StaticOn();
  renderer->AddActor(lodActor.GetPointer());

  vtkNew<vtkPolyDataMapper> wmapper;
  wmapper->SetInputData(pd);
  wmapper->ScalarVisibilityOff();

  vtkNew<vtkActor> wactor;
  wactor->SetMapper(wmapper.GetPointer());
  wactor->GetProperty()->SetColor(0, 0, 0);
  wactor->GetProperty()->SetRepresentationToWireframe();
  //wactor->SetScale(1,1,2);
  //renderer->AddActor(wactor.GetPointer());

  vtkNew<vtkPolyDataMapper> slmapper;
  slmapper->SetInputConnection(0, slsurface->GetOutputPort(0));
  slmapper->SetScalarRange(dataRange);

  vtkNew<vtkActor> slactor;
  slactor->SetMapper(slmapper.GetPointer());
  slactor->SetPosition(0, 0, 1500);
  renderer->AddActor(slactor.GetPointer());

  vtkNew<vtkPolyDataMapper> slwmapper;
  slwmapper->SetInputConnection(0, slsurface->GetOutputPort(0));
  slwmapper->ScalarVisibilityOff();

  vtkNew<vtkActor> slwactor;
  slwactor->SetMapper(slwmapper.GetPointer());
  slwactor->SetPosition(0, 0, 1500);
  slwactor->GetProperty()->SetColor(0, 0, 0);
  slwactor->GetProperty()->SetRepresentationToWireframe();
  renderer->AddActor(slwactor.GetPointer());

  vtkNew<vtkAxesActor> axes;
  axes->SetTotalLength(5000, 5000, 5000);
  //renderer->AddActor(axes.GetPointer());

  // Standard testing code.
  //renWin->SetMultiSamples(0);
  renderer->SetBackground(1., 1., 1.);
  renWin->AddRenderer(renderer.GetPointer());
  renWin->SetSize(900, 900);

  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin.GetPointer());
  iren->SetDesiredUpdateRate(25.);
  iren->SetStillUpdateRate(0.001);
  renWin->SetDesiredUpdateRate(25.);
  renWin->Render();

  vtkNew<MytInteractorStyle> style;
  iren->SetInteractorStyle(style.GetPointer());
  iren->SetRenderWindow(renWin.GetPointer());
  iren->Initialize();
  iren->Start();
  //int retVal = vtkRegressionTestImage(renWin.GetPointer());
 // if ( retVal == vtkRegressionTester::DO_INTERACTOR)
    {
   // iren->Start();
    }

  return EXIT_SUCCESS; //!retVal;
}
