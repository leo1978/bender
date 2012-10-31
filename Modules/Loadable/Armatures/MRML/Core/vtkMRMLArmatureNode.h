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

#ifndef __vtkMRMLArmatureNode_h
#define __vtkMRMLArmatureNode_h

// Slicer includes
#include "vtkMRMLAnnotationHierarchyNode.h"

// Armatures includes
#include "vtkBenderArmaturesModuleMRMLCoreExport.h"
class vtkMRMLBoneNode;

/// \ingroup Bender_MRML
/// \brief Root of a tree of bones
///
/// \sa vtkMRMLBoneNode
class VTK_BENDER_ARMATURES_MRML_CORE_EXPORT vtkMRMLArmatureNode
  : public vtkMRMLAnnotationHierarchyNode
{
public:
  //--------------------------------------------------------------------------
  // VTK methods
  //--------------------------------------------------------------------------

  static vtkMRMLArmatureNode *New();
  vtkTypeMacro(vtkMRMLArmatureNode,vtkMRMLAnnotationHierarchyNode);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------

  /// Instantiate a bone node.
  virtual vtkMRMLNode* CreateNodeInstance();

  /// Get node XML tag name (like Volume, Model).
  virtual const char* GetNodeTagName() {return "Armature";};

  /// Read node attributes from XML file.
  virtual void ReadXMLAttributes( const char** atts);

  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  /// Copy the node's attributes to this object.
  virtual void Copy(vtkMRMLNode *node);

  /// Update references from scene.
  virtual void UpdateScene(vtkMRMLScene *scene);

  /// Alternative method to propagate events generated by observed nodes.
  virtual void ProcessMRMLEvents(vtkObject* caller,
                                 unsigned long event,
                                 void* callData);

  //--------------------------------------------------------------------------
  // Annotation methods
  //--------------------------------------------------------------------------

  virtual const char* GetIcon();

  //--------------------------------------------------------------------------
  // Bone methods
  //--------------------------------------------------------------------------
  /// Search and populate the \a bones collection with all the bones found
  /// in the armature.
  /// \sa GetBoneParent()
  inline void GetAllBones(vtkCollection* bones);

  ///.Search and return the parent of the \a boneNode bone. Return 0 if the
  /// bone doesn't have a parent (top-level) or doesn't belong to the
  /// armature.
  /// \sa GetAllBones()
  vtkMRMLBoneNode* GetParentBone(vtkMRMLBoneNode* boneNode);
protected:
  vtkMRMLArmatureNode();
  ~vtkMRMLArmatureNode();

  vtkMRMLArmatureNode(const vtkMRMLArmatureNode&); /// not implemented
  void operator=(const vtkMRMLArmatureNode&); /// not implemented
};

//----------------------------------------------------------------------------
void vtkMRMLArmatureNode::GetAllBones(vtkCollection* bones)
{
  this->GetAllChildren(bones);
}

#endif
