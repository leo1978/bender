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

#ifndef __vtkMRMLBoneNode_h
#define __vtkMRMLBoneNode_h

// Slicer includes
#include "vtkMRMLAnnotationNode.h"

// Armatures includes
#include "vtkBenderArmaturesModuleMRMLCoreExport.h"
class vtkMRMLBoneDisplayNode;
class vtkBoneWidget;

/// \ingroup Bender_MRML
/// \brief Annotation to design and edit a bone.
///
/// In comparison with annotation nodes, vtkMRMLBoneNode only supports
/// vtkMRMLBoneDisplayNode for display nodes.
/// \sa vtkMRMLBoneDisplayNode, vtkMRMLArmatureNode
class VTK_BENDER_ARMATURES_MRML_CORE_EXPORT vtkMRMLBoneNode
  : public vtkMRMLAnnotationNode
{
public:
  //--------------------------------------------------------------------------
  // VTK methods
  //--------------------------------------------------------------------------

  static vtkMRMLBoneNode *New();
  vtkTypeMacro(vtkMRMLBoneNode,vtkMRMLAnnotationNode);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------

  /// Instantiate a bone node.
  virtual vtkMRMLNode* CreateNodeInstance();

  /// Get node XML tag name (like Volume, Model).
  virtual const char* GetNodeTagName() {return "Bone";};

  /// Read node attributes from XML file.
  virtual void ReadXMLAttributes( const char** atts);

  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  /// Copy the node's attributes to this object.
  virtual void Copy(vtkMRMLNode *node) {Superclass::Copy(node);}

  /// Update references from scene.
  virtual void UpdateScene(vtkMRMLScene *scene);

  /// Alternative method to propagate events generated by observed nodes.
  virtual void ProcessMRMLEvents(vtkObject* caller,
                                 unsigned long event,
                                 void* callData);

  //--------------------------------------------------------------------------
  // Annotation methods
  //--------------------------------------------------------------------------

  virtual void Initialize(vtkMRMLScene* mrmlScene);

  //--------------------------------------------------------------------------
  // Bone methods
  //--------------------------------------------------------------------------

  /// Utility function that returns the associated display node as a bone
  /// display node.
  /// \sa GetDisplayNode()
  vtkMRMLBoneDisplayNode* GetBoneDisplayNode();
  /// Create a default display node if not already present.
  /// \sa CreateDefaultStorageNode()
  void CreateBoneDisplayNode();

  /// Set the bone head position in world coordinates.
  /// \sa GetWorldHeadRest(), SetWorldTailRest()
  void SetWorldHeadRest(double headPoint[3]);
  /// Get the bone head position in world coordinates.
  /// \sa SetWorldHeadRest(), GetWorldTailRest()
  double* GetWorldHeadRest();

  /// Set the bone tail position in world coordinates.
  /// \sa GetWorldTailRest(), SetWorldHeadRest()
  void SetWorldTailRest(double tailPoint[3]);
  /// Get the bone tail position in world coordinates.
  /// \sa GetWorldTailRest(), GetWorldHeadRest()
  double* GetWorldTailRest();

protected:
  vtkMRMLBoneNode();
  ~vtkMRMLBoneNode();

  vtkMRMLBoneNode(const vtkMRMLBoneNode&); /// not implemented
  void operator=(const vtkMRMLBoneNode&); /// not implemented

  vtkBoneWidget* BoneProperties;
};

#endif