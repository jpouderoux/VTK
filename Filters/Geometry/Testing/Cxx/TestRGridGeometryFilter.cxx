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
#include "vtkCellData.h"
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

  pts->SetNumberOfPoints(si * sj * sk * 8);
  for (int k = 0, ijk = 0; k < sk; k++)
    {
    for (int j = 0; j < sj; j++)
      {
      for (int i = 0; i < si; i++, ijk++)
        {
        double z = i < 30 ? k : k + 0.7;

        pts->SetPoint(ijk * 8 + 0, i, j, z);
        pts->SetPoint(ijk * 8 + 1, i+1, j, z);
        pts->SetPoint(ijk * 8 + 2, i+1, j+1, z);
        pts->SetPoint(ijk * 8 + 3, i, j+1, z);

        pts->SetPoint(ijk * 8 + 4, i, j, z+1);
        pts->SetPoint(ijk * 8 + 5, i+1, j, z+1);
        pts->SetPoint(ijk * 8 + 6, i+1, j+1, z+1);
        pts->SetPoint(ijk * 8 + 7, i, j+1, z+1);
        }
      }
    }

  vtkNew<vtkRGrid> grid;
  grid->SetDimensions(si, sj, sk);
  grid->SetPoints(pts.GetPointer());

  cout << grid->GetNumberOfCells() << " cells "
    << grid->GetNumberOfPoints() << " pts "<< endl;

  grid->PrintSelf(cout, vtkIndent());
  //vtkCell* c = grid->GetCell(1);
  //c->PrintSelf(cout, vtkIndent());
  cout << "---------------" << endl;

  vtkNew<vtkDataSetSurfaceFilter> surface;
  surface->SetInputData(grid.GetPointer());
  surface->Update();

  vtkPolyData *pd = surface->GetOutput();

  pd->PrintSelf(cout, vtkIndent());
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

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());
  renderer->AddActor(actor.GetPointer());

  // Standard testing code.
  renWin->SetMultiSamples(0);
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
