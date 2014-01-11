/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestRGrid.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkNew.h"
#include "vtkRGrid.h"
#include "vtkPoints.h"
#include "vtkCell.h"

int TestRGrid(int,char *[])
{
  vtkNew<vtkRGrid> grid;

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
  grid->SetDimensions(si, sj, sk);
  grid->SetPoints(pts.GetPointer());

  vtkCell* c = grid->GetCell(1);
  c->PrintSelf(std::cout, vtkIndent());


  return EXIT_SUCCESS;
}
