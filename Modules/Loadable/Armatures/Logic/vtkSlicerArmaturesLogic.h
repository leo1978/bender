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

// .NAME vtkSlicerArmaturesLogic - slicer logic class for armature manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing properties of armatures

#ifndef __vtkSlicerArmaturesLogic_h
#define __vtkSlicerArmaturesLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"
#include "vtkSlicerArmaturesModuleLogicExport.h"
class vtkSlicerModelsLogic;

// MRML includes
class vtkMRMLModelNode;

// VTK includes
class vtkPolyData;
class vtkXMLDataElement;
class vtkXMLDataParser;

// STD includes
#include <cstdlib>

/// \ingroup Slicer_QtModules_Armatures
class VTK_SLICER_ARMATURES_MODULE_LOGIC_EXPORT vtkSlicerArmaturesLogic
  : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerArmaturesLogic *New();
  vtkTypeMacro(vtkSlicerArmaturesLogic,vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkMRMLModelNode* AddArmatureFile(const char* fileName, bool applyPose = true);

  virtual void SetModelsLogic(vtkSlicerModelsLogic* modelsLogic);
  vtkGetObjectMacro(ModelsLogic, vtkSlicerModelsLogic);

  vtkGetMacro(LPS2RAS, double);
  vtkSetMacro(LPS2RAS, double);

  static void ComputeTransform(double start[3], double end[3], double mat[3][3]);
  static double ComputeAngle(double v1[3], double v2[3]);
  static void ComputeAxisAngleMatrix(double axis[3], double angle, double mat[3][3]);

protected:
  vtkSlicerArmaturesLogic();
  virtual ~vtkSlicerArmaturesLogic();

  void ReadArmature(vtkXMLDataParser* parser, vtkPolyData* polyData, bool applyPose, bool applyVertexGroup);
  void ReadBone(vtkXMLDataElement* boneElement, vtkPolyData* polyData, bool applyPose, bool applyVertexGroup,
                const double origin[3], const double parentMat[3][3], double parentLength = 0.);
  void ReadPose(vtkXMLDataElement* poseElement, double orientation[4], bool invert = true);

  vtkSlicerModelsLogic* ModelsLogic;
  double LPS2RAS;

private:
  vtkSlicerArmaturesLogic(const vtkSlicerArmaturesLogic&); // Not implemented
  void operator=(const vtkSlicerArmaturesLogic&);               // Not implemented
};

#endif

