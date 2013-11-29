/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPPCAStatistics.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2011 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
  -------------------------------------------------------------------------*/
#include "vtkToolkits.h"

#include "vtkPPCAStatistics.h"

#include "vtkAbstractArray.h"
#include "vtkCommunicator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiProcessController.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPOrderStatistics.h"
#include "vtkTable.h"
#include "vtkPMultiCorrelativeStatistics.h"
#include "vtkVariant.h"

#include <map>

vtkStandardNewMacro(vtkPPCAStatistics);
vtkCxxSetObjectMacro(vtkPPCAStatistics, Controller, vtkMultiProcessController);
//-----------------------------------------------------------------------------
vtkPPCAStatistics::vtkPPCAStatistics()
{
  this->Controller = 0;
  this->SetController( vtkMultiProcessController::GetGlobalController() );
}

//-----------------------------------------------------------------------------
vtkPPCAStatistics::~vtkPPCAStatistics()
{
  this->SetController( 0 );
}

//-----------------------------------------------------------------------------
void vtkPPCAStatistics::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Controller: " << this->Controller << endl;
}

// ----------------------------------------------------------------------
void vtkPPCAStatistics::Learn( vtkTable* inData,
                               vtkTable* inParameters,
                               vtkMultiBlockDataSet* outMeta )
{
  if ( ! outMeta )
    {
    return;
    }
    
  // First calculate correlative statistics on local data set
  this->Superclass::Learn( inData, inParameters, outMeta );

  // Get a hold of the (sparse) covariance matrix
  vtkTable* sparseCov = vtkTable::SafeDownCast( outMeta->GetBlock( 0 ) );
  if ( ! sparseCov )
    {
    return;
    }

  if ( !this->MedianAbosluteVariance )
    {
    vtkPMultiCorrelativeStatistics::GatherStatistics( this->Controller, sparseCov );
    }
}

// ----------------------------------------------------------------------
void vtkPPCAStatistics::Test( vtkTable* inData,
                              vtkMultiBlockDataSet* inMeta,
                              vtkTable* outMeta )
{
  if ( this->Controller->GetNumberOfProcesses() > 1 )
    {
    vtkWarningMacro( "Parallel PCA: Hypothesis testing not implemented for more than 1 process." );
    return;
    }

  this->Superclass::Test( inData, inMeta, outMeta );
}

void vtkPPCAStatistics::ComputeMedian(vtkTable* inData, vtkTable* outData)
{
  vtkNew<vtkPOrderStatistics> orderStats;
  vtkNew<vtkTable> inOrderStats;
  orderStats->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, inOrderStats.GetPointer());
  for (vtkIdType i = 0; i < inData->GetNumberOfColumns(); ++ i )
    {
    inOrderStats->AddColumn(inData->GetColumn(i));
    orderStats->AddColumn(inData->GetColumn(i)->GetName());
    }
  orderStats->SetNumberOfIntervals(2);
  orderStats->SetLearnOption(true);
  orderStats->SetDeriveOption(true);
  orderStats->SetTestOption(false);
  orderStats->SetAssessOption(false);
  orderStats->Update();
  // Gets the Median
  vtkMultiBlockDataSet *outputOrderStats =
    vtkMultiBlockDataSet::SafeDownCast(
    orderStats->GetOutputDataObject(vtkStatisticsAlgorithm::OUTPUT_MODEL));
  outData->DeepCopy(vtkTable::SafeDownCast(
    outputOrderStats->GetBlock(outputOrderStats->GetNumberOfBlocks() - 1)));

   return ;
}
