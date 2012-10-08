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

#ifndef __vtkArmatureWidget_h
#define __vtkArmatureWidget_h

// .NAME vtkArmatureWidget - Widget for holding bone widgets together
// .SECTION Description
// .SECTION See Also
// vtkArmatureRepresentation, vtkBoneWidget
// Bender includes
#include "vtkBenderWidgetsExport.h"

// VTK includes
#include <vtkAbstractWidget.h>

class vtkArmatureRepresentation;
class vtkArmatureWidgetCallback;
class vtkBoneWidget;

class VTK_BENDER_WIDGETS_EXPORT vtkArmatureWidget : public vtkAbstractWidget
{
public:
  // Description:
  // Instantiate this class.
  static vtkArmatureWidget *New();

  // Description:
  // Standard methods for a VTK class.
  vtkTypeMacro(vtkArmatureWidget, vtkAbstractWidget);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The method for activiating and deactiviating this widget. This method
  // must be overridden because it is a composite widget and does more than
  // its superclasses' vtkAbstractWidget::SetEnabled() method.
  virtual void SetEnabled(int);

  // Description:
  void SetRepresentation(vtkArmatureRepresentation* r);

  // Description:
  // Return the representation as a vtkArmatureRepresentation.
  vtkArmatureRepresentation* GetArmatureRepresentation();

  // Description:
  // Create the default widget representation (vtkArmatureRepresentation) if no one is set.
  virtual void CreateDefaultRepresentation();

  // Description:
  // Methods to change the whether the widget responds to interaction.
  // Overridden to pass the state to bone widgets.
  virtual void SetProcessEvents(int process);

protected:
  vtkArmatureWidget();
  ~vtkArmatureWidget();

//BTX
  class vtkBoneWidgetCallback;
//ETX
  vtkBoneWidgetCallback* BoneWidgetCallback;

  vtkBoneWidget* RootBoneWidget;

private:
  vtkArmatureWidget(const vtkArmatureWidget&);  //Not implemented
  void operator=(const vtkArmatureWidget&);  //Not implemented
};

#endif
