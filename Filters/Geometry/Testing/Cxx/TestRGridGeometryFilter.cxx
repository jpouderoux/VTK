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
#include <cassert>
#include "vtkLookupTable.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGridGeometryFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkGeometryFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkTimerLog.h"
#include "vtkCamera.h"
#include "vtkPointData.h"
#include "vtkProperty.h"
#include "vtkRGrid.h"
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
  const int si1 = si+1;
  const int sj1 = sj+2;
  const int sk1 = sk+2;

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
          z = k + 0.7;
          if ( i > 50)
            {
            x = x-1;
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
  scalars->SetNumberOfValues(si*sj*sk);

  vtkNew<vtkCellArray> cells;
  cells->Allocate(9 * si * sj * sk);

  for (int k = 0, ijk = 0; k < sk; k++)
    {
    for (int j = 0; j < sj; j++)
      {
      for (int i = 0; i < si; i++, ijk++)
        {
        cells->InsertNextCell(8);
        cells->InsertCellPoint(i   + j     * si1 + k * sj1 * si1);
        cells->InsertCellPoint(i+1 + j     * si1 + k * sj1 * si1);
        cells->InsertCellPoint(i+1 + (j+1) * si1 + k * sj1 * si1);
        cells->InsertCellPoint(i   + (j+1) * si1 + k * sj1 * si1);

        cells->InsertCellPoint(i   + j     * si1 + (k+1) * sj1 * si1);
        cells->InsertCellPoint(i+1 + j     * si1 + (k+1) * sj1 * si1);
        cells->InsertCellPoint(i+1 + (j+1) * si1 + (k+1) * sj1 * si1);
        cells->InsertCellPoint(i   + (j+1) * si1 + (k+1) * sj1 * si1);
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
  vtkNew<vtkIdList> ids;
  grid->GetPointCells(62+(61*51), ids.GetPointer());
  //grid->PrintSelf(cout, vtkIndent());
  for (int i=0; i < ids->GetNumberOfIds(); i++)
    cout << i << ": " << ids->GetId(i) << endl;
  ids->PrintSelf(cout, vtkIndent());
  vtkCell* c = grid->GetCell(0);
  c->PrintSelf(cout, vtkIndent());
  c = grid->GetCell(1);
  c->PrintSelf(cout, vtkIndent());
  cout << "---------------" << endl;

  vtkNew<vtkDataSetSurfaceFilter> surface;
  surface->SetInputData(grid.GetPointer());
  surface->Update();

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
  wactor->GetProperty()->SetColor(0,0,0);
  wactor->GetProperty()->SetRepresentationToWireframe();
  renderer->AddActor(wactor.GetPointer());

  // Standard testing code.
  //renWin->SetMultiSamples(0);
  renderer->SetBackground(0.5, 0.5, 0.5);
  renWin->SetSize(300, 300);
  renWin->Render();
  int retVal = vtkRegressionTestImage(renWin.GetPointer());
 // if ( retVal == vtkRegressionTester::DO_INTERACTOR)
    {
    iren->Start();
    }

  return !retVal;
}
