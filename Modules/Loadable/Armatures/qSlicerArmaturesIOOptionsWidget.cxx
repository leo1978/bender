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

/// Qt includes
#include <QFileInfo>

// CTK includes
#include <ctkFlowLayout.h>
#include <ctkUtils.h>

/// Armatures includes
#include "qSlicerIOOptions_p.h"
#include "qSlicerArmaturesIOOptionsWidget.h"
#include "ui_qSlicerArmaturesIOOptionsWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_Armatures
class qSlicerArmaturesIOOptionsWidgetPrivate
  : public qSlicerIOOptionsPrivate
  , public Ui_qSlicerArmaturesIOOptionsWidget
{
public:
};

//-----------------------------------------------------------------------------
qSlicerArmaturesIOOptionsWidget::qSlicerArmaturesIOOptionsWidget(QWidget* parentWidget)
  : qSlicerIOOptionsWidget(new qSlicerArmaturesIOOptionsWidgetPrivate, parentWidget)
{
  Q_D(qSlicerArmaturesIOOptionsWidget);
  d->setupUi(this);

  connect(d->PoseCheckBox, SIGNAL(toggled(bool)),
          this, SLOT(updateProperties()));

  // Apply pose by default
  d->PoseCheckBox->setChecked(true);
}

//-----------------------------------------------------------------------------
qSlicerArmaturesIOOptionsWidget::~qSlicerArmaturesIOOptionsWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerArmaturesIOOptionsWidget::updateProperties()
{
  Q_D(qSlicerArmaturesIOOptionsWidget);
  d->Properties["applyPose"] = d->PoseCheckBox->isChecked();
}

