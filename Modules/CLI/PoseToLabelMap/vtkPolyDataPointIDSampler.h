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

#ifndef __vtkPolyDataPointIDSampler_h
#define __vtkPolyDataPointIDSampler_h

#include "vtkPolyDataPointSampler.h"

class vtkPolyDataPointIDSampler
  : public vtkPolyDataPointSampler
{
public:
  // Description:
  // Instantiate this class.
  static vtkPolyDataPointIDSampler *New();

  // Description:
  // Standard macros for type information and printing.
  vtkTypeMacro(vtkPolyDataPointIDSampler,vtkPolyDataPointSampler);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkPolyDataPointIDSampler();
  ~vtkPolyDataPointIDSampler();

  virtual void InitializeOutput(vtkPolyData* output, vtkPolyData* input);
  virtual void InsertNextPoint(vtkPolyData *output, double x[3],
    vtkPoints *inPts, vtkIdType npts, vtkIdType *pts,
    double s, double t = 0.);

private:
  vtkPolyDataPointIDSampler(const vtkPolyDataPointIDSampler&);  // Not implemented.
  void operator=(const vtkPolyDataPointIDSampler&);  // Not implemented.
};

#endif
