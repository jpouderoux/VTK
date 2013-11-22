#include "vtkComputeQuartiles.h"
#include "vtkStatisticsAlgorithm.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkNew.h"
#include "vtkTable.h"

//----------------------------------------------------------------------------
int TestComputeQuartiles(int , char * [])
{
  vtkNew<vtkDoubleArray> arrFirstVariable;
  arrFirstVariable->SetName("Math");
  
  vtkNew<vtkDoubleArray> arrSecondVariable;
  arrSecondVariable->SetName("French");
    
  // Create a table with two columns in it...
  vtkNew<vtkTable> table;
  table->AddColumn(arrFirstVariable.GetPointer());
  table->AddColumn(arrSecondVariable.GetPointer());
  
  const int numNotes = 20;
  table->SetNumberOfRows(numNotes);

  const double MathValue[] =
    {
    18, 20, 20, 16,
    12, 14, 16, 14,
    14, 13, 16, 18,
    6, 10, 16, 14,
    4, 16, 16, 14
    };

  const double FrenchValue[] =
    {
    14, 12, 14, 16,
    12, 14, 16, 4,
    4, 10, 6, 20,
    14, 16, 14, 14,
    12, 2, 14, 8
    };

  for (int i = 0; i < numNotes; ++i)
    {
    table->SetValue(i, 0, MathValue[i]);
    table->SetValue(i, 1, FrenchValue[i]);
    }

  // Run Compute Quantiles
  
  vtkNew<vtkComputeQuartiles> quartiles;

  // First verify that absence of input does not cause trouble
  cout << "## Verifying that absence of input does not cause trouble... ";
  quartiles->Update();
  cout << "done.\n";

  quartiles->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, table.GetPointer());
  quartiles->Update();

  vtkTable* quartileTables = quartiles->GetOutput();
  quartileTables->Dump();

  return EXIT_SUCCESS;
}
