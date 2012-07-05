/*=========================================================================

  Program: Bender

  Copyright (c) Kitware Inc.

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=========================================================================*/

#include "vtkPolyDataPointIDSampler.h"

// VTK includes
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

vtkStandardNewMacro(vtkPolyDataPointIDSampler);

//------------------------------------------------------------------------
vtkPolyDataPointIDSampler::vtkPolyDataPointIDSampler()
{
}

//------------------------------------------------------------------------
vtkPolyDataPointIDSampler::~vtkPolyDataPointIDSampler()
{
}

//------------------------------------------------------------------------
void vtkPolyDataPointIDSampler
::InitializeOutput(vtkPolyData *output, vtkPolyData* input)
{
  this->Superclass::InitializeOutput(output, input);
  vtkNew<vtkIdTypeArray> ids;
  ids->SetName("pointIndexes");
  output->GetPointData()->SetScalars(ids.GetPointer());
}

//------------------------------------------------------------------------
void vtkPolyDataPointIDSampler
::InsertNextPoint(vtkPolyData *output, double x[3], vtkPoints *inPts,
                  vtkIdType npts, vtkIdType *pts,
                  double s, double t)
{
  output->GetPoints()->InsertNextPoint(x);
  vtkIdTypeArray* ids = vtkIdTypeArray::SafeDownCast(
    output->GetPointData()->GetScalars("pointIndexes"));
  vtkIdType closestpointId = s >  0.5 ? (t > 0.5 ? pts[2] : pts[1]) :
                             (t > 0.5 ? pts[ (npts > 3 ? 3 : 2)] : pts[0]);
  ids->InsertNextValue(closestpointId);
}

//------------------------------------------------------------------------
void vtkPolyDataPointIDSampler::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
