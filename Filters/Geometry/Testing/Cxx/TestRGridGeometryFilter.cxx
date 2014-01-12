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
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkLookupTable.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGridGeometryFilter.h"
#include "vtkDataSetSurfaceFilter.h"
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

int TestRGridGeometryFilter(int argc, char* argv[])
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
          z = k + 1.0;
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
        }
      }
    }

  vtkNew<vtkRGrid> grid;
  grid->SetDimensions(si, sj, sk);
  grid->SetPoints(pts.GetPointer());
  grid->SetCells(cells.GetPointer());
  grid->GetCellData()->AddArray(scalars.GetPointer());
  grid->GetCellData()->SetActiveScalars("scalars");

  cout << grid->GetNumberOfCells() << " cells "
    << grid->GetNumberOfPoints() << " pts "<< endl;

 vtkNew<vtkRGridSlice> slice;
  slice->SetInputData(grid.GetPointer());
  slice->Update();

  slice->GetOutput()->PrintSelf(cout, vtkIndent());
  vtkRGrid* sl = slice->GetOutput();

  vtkNew<vtkDataSetSurfaceFilter> surface;
  surface->SetInputData(grid.GetPointer());
  surface->Update();

  vtkNew<vtkDataSetSurfaceFilter> slsurface;
  slsurface->SetInputData(sl);
  slsurface->Update();

  vtkPolyData *pd = surface->GetOutput();

  //pd->PrintSelf(cout, vtkIndent());
  cout << pd->GetNumberOfCells() << " cells "
    << pd->GetNumberOfPoints() << " pts "<< endl;

  vtkNew<vtkOutlineCornerFilter> outline;
  outline->SetInputData(grid.GetPointer());

  // Standard rendering classes
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(renderer.GetPointer());
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin.GetPointer());


  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(0, surface->GetOutputPort(0));
  mapper->SetScalarRange(0, sk);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());
  renderer->AddActor(actor.GetPointer());

  vtkNew<vtkPolyDataMapper> wmapper;
  wmapper->SetInputConnection(0, surface->GetOutputPort(0));
  wmapper->ScalarVisibilityOff();

  vtkNew<vtkActor> wactor;
  wactor->SetMapper(wmapper.GetPointer());
  wactor->GetProperty()->SetColor(0, 0, 0);
  wactor->GetProperty()->SetRepresentationToWireframe();
  renderer->AddActor(wactor.GetPointer());


  vtkNew<vtkPolyDataMapper> slmapper;
  slmapper->SetInputConnection(0, slsurface->GetOutputPort(0));
  slmapper->SetScalarRange(0, sk);

  vtkNew<vtkActor> slactor;
  slactor->SetMapper(slmapper.GetPointer());
  slactor->SetPosition(0,0,25);
  renderer->AddActor(slactor.GetPointer());

  vtkNew<vtkPolyDataMapper> slwmapper;
  slwmapper->SetInputConnection(0, slsurface->GetOutputPort(0));
  slwmapper->ScalarVisibilityOff();

  vtkNew<vtkActor> slwactor;
  slwactor->SetMapper(slwmapper.GetPointer());
  slwactor->SetPosition(0,0,25);
  slwactor->GetProperty()->SetColor(0, 0, 0);
  slwactor->GetProperty()->SetRepresentationToWireframe();
  renderer->AddActor(slwactor.GetPointer());


  // Standard testing code.
  //renWin->SetMultiSamples(0);
  renderer->SetBackground(0.5, 0.5, 0.5);
  renWin->SetSize(900, 900);
  renWin->Render();
  int retVal = vtkRegressionTestImage(renWin.GetPointer());
 // if ( retVal == vtkRegressionTester::DO_INTERACTOR)
    {
    iren->Start();
    }

  return !retVal;
}
