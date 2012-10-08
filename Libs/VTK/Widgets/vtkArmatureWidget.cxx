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

// Bender includes
#include "vtkArmatureRepresentation.h"
#include "vtkArmatureWidget.h"
#include "vtkBoneWidget.h"

// VTK includes
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkWidgetCallbackMapper.h>
#include <vtkWidgetEvent.h>

vtkStandardNewMacro(vtkArmatureWidget);

class vtkArmatureWidget::vtkBoneWidgetCallback : public vtkCommand
{
public:
  static vtkBoneWidgetCallback *New() {return new vtkBoneWidgetCallback;}
  vtkBoneWidgetCallback() { this->ArmatureWidget = 0; }
  virtual void Execute(vtkObject* caller, unsigned long eventId, void*)
    {
    switch (eventId)
      {
      case vtkCommand::StartInteractionEvent:
        {
        //this->ArmatureWidget->StartBoneInteraction();
        break;
        }
      case vtkCommand::EndInteractionEvent:
        {
        //this->ArmatureWidget->EndBoneInteraction();
        break;
        }
      }
    }
  vtkArmatureWidget *ArmatureWidget;
};

//----------------------------------------------------------------------
vtkArmatureWidget::vtkArmatureWidget()
{
  this->BoneWidgetCallback = vtkArmatureWidget::vtkBoneWidgetCallback::New();
  this->BoneWidgetCallback->ArmatureWidget = this;

  this->RootBoneWidget = vtkBoneWidget::New();
  this->RootBoneWidget->SetParent(this);
  this->RootBoneWidget->ManagesCursorOff();

  this->RootBoneWidget->AddObserver(vtkCommand::StartInteractionEvent, this->BoneWidgetCallback,
                                    this->Priority);
  this->RootBoneWidget->AddObserver(vtkCommand::EndInteractionEvent, this->BoneWidgetCallback,
                                    this->Priority);
}

//----------------------------------------------------------------------
vtkArmatureWidget::~vtkArmatureWidget()
{
  this->RootBoneWidget->RemoveObserver(this->BoneWidgetCallback);
  this->RootBoneWidget->Delete();
  this->BoneWidgetCallback->Delete();
}

//----------------------------------------------------------------------
void vtkArmatureWidget::CreateDefaultRepresentation()
{
  vtkArmatureRepresentation* armatureRepresentation =
    vtkArmatureRepresentation::SafeDownCast(this->WidgetRep);
  if (!armatureRepresentation)
    {
    armatureRepresentation = vtkArmatureRepresentation::New();
    this->WidgetRep = armatureRepresentation;
    }
  //armatureRepresentation->InstantiateHandleRepresentation();
}

//----------------------------------------------------------------------
void vtkArmatureWidget::SetEnabled(int enabling)
{
  if (enabling)
    {
    this->RootBoneWidget->SetInteractor(this->Interactor);
    this->RootBoneWidget->GetRepresentation()->SetRenderer(
      this->CurrentRenderer);
    }
  if (this->WidgetRep)
    {
    this->WidgetRep->SetVisibility(enabling);
    }
  this->RootBoneWidget->SetEnabled(enabling);
  this->Superclass::SetEnabled(enabling);
}

//----------------------------------------------------------------------
void vtkArmatureWidget::SetRepresentation(vtkArmatureRepresentation* r)
{
  this->Superclass::SetWidgetRepresentation(r);
}

//----------------------------------------------------------------------
vtkArmatureRepresentation* vtkArmatureWidget::GetArmatureRepresentation()
{
  return vtkArmatureRepresentation::SafeDownCast(this->WidgetRep);
}

//----------------------------------------------------------------------
void vtkArmatureWidget::SetProcessEvents(int pe)
{
  this->Superclass::SetProcessEvents(pe);

  this->RootBoneWidget->SetProcessEvents(pe);
}

//----------------------------------------------------------------------
void vtkArmatureWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Armature Widget " << this << "\n";
  os << indent << "Root Bone: "<< this->RootBoneWidget<< "\n";
}
