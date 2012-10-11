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

#ifndef __qSlicerArmaturesIOOptionsWidget_h
#define __qSlicerArmaturesIOOptionsWidget_h

// SlicerQt includes
#include "qSlicerIOOptionsWidget.h"

// Armatures includes
#include "qSlicerArmaturesModuleExport.h"
class qSlicerArmaturesIOOptionsWidgetPrivate;

/// \ingroup Slicer_QtModules_Armatures
class Q_SLICER_QTMODULES_ARMATURES_EXPORT qSlicerArmaturesIOOptionsWidget :
  public qSlicerIOOptionsWidget
{
  Q_OBJECT
public:
  qSlicerArmaturesIOOptionsWidget(QWidget *parent=0);
  virtual ~qSlicerArmaturesIOOptionsWidget();

protected slots:
  void updateProperties();

private:
  Q_DECLARE_PRIVATE_D(qGetPtrHelper(qSlicerIOOptions::d_ptr), qSlicerArmaturesIOOptionsWidget);
  Q_DISABLE_COPY(qSlicerArmaturesIOOptionsWidget);
};

#endif
