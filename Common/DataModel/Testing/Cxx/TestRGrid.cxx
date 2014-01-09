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


int TestRGrid(int,char *[])
{
  vtkNew<vtkRGrid> grid;


  return EXIT_SUCCESS;
}
