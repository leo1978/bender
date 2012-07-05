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

// ModelToLabelMap includes
#include "PoseToLabelMapCLP.h"
#include "vtkPolyDataPointIDSampler.h"

// ITK includes
#ifdef ITKV3_COMPATIBILITY
#include <itkAnalyzeImageIOFactory.h>
#endif
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkPluginUtilities.h>

// VTK includes
#include <vtkDebugLeaks.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>

//-----------------------------------------------------------------------------
template<class T>
struct InputImage{
  typedef itk::Image<T, 3> Type;
};
typedef itk::Image<unsigned int, 3> VoxelizedModelImageType;

template <class T>
int DoIt( int argc, char * argv[] );
template <class T>
VoxelizedModelImageType::Pointer VoxelizeModel(vtkPolyData* model,
                                               typename InputImage<T>::Type::Pointer inputImage,
                                               double samplingDistance);
vtkPolyData* ReadPolyData(const std::string& fileName);
template<class T>
void WriteImage(typename InputImage<T>::Type::Pointer image,
                const std::string& fileName,
                ModuleProcessInformation* processInformation);

//-----------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;
#ifdef ITKV3_COMPATIBILITY
  itk::ObjectFactoryBase::RegisterFactory( itk::AnalyzeImageIOFactory::New() );
#endif
  try
    {
    itk::GetImageType(InputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types

    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
      case itk::ImageIOBase::CHAR:
        return DoIt<char>( argc, argv );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::SHORT:
        return DoIt<short>( argc, argv );
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::INT:
        return DoIt<int>( argc, argv );
        break;
      case itk::ImageIOBase::ULONG:
      case itk::ImageIOBase::LONG:
        return DoIt<long>( argc, argv );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt<float>( argc, argv );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt<double>( argc, argv );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
template <class T>
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;
  vtkDebugLeaks::SetExitError(true);

  typedef T InputPixelType;

  typedef itk::ImageFileReader<typename InputImage<T>::Type> ReaderType;
  typedef itk::ImageFileWriter<VoxelizedModelImageType> WriterType;

  // Read the input volume
  typename ReaderType::Pointer reader = ReaderType::New();
  itk::PluginFilterWatcher watchReader(reader, "Read Input Volume",
                                       CLPProcessInformation);
  reader->SetFileName( InputVolume.c_str() );
  reader->Update();

  vtkPolyData* polyData = ReadPolyData(InputModel);
  if (polyData == 0)
    {
    return EXIT_FAILURE;
    }

// read the poly data
  vtkSmartPointer<vtkPolyData> model;
  model.TakeReference(polyData);
  
  // output label map
  VoxelizedModelImageType::Pointer label =
    VoxelizeModel<T>(model, reader->GetOutput(), samplingDistance);

  if (debug)
    {
    WriteImage<VoxelizedModelImageType::ValueType>(label, "e:\\voxelizedMap.mha", CLPProcessInformation);
    }

  // Distance map
  typedef itk::Image<unsigned short, 3> DistanceMapImageType;
  typedef itk::DanielssonDistanceMapImageFilter<VoxelizedModelImageType, DistanceMapImageType > DistanceMapFilterType;
  DistanceMapFilterType::Pointer filter = DistanceMapFilterType::New();
  filter->SetInput( label );
  filter->Update();
  DistanceMapImageType::Pointer distanceMap = filter->GetOutput();

  if (debug)
    {
    WriteImage<DistanceMapImageType::ValueType>(distanceMap,
      "e:\\distanceMap.mha",
      CLPProcessInformation);
    }

  //

  WriteImage<DistanceMapImageType::ValueType>(distanceMap,
    OutputVolume.c_str(),
    CLPProcessInformation);

  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
template <class T>
VoxelizedModelImageType::Pointer VoxelizeModel(vtkPolyData* model,
                                               typename InputImage<T>::Type::Pointer inputImage,
                                               double samplingDistance)
{
  VoxelizedModelImageType::Pointer label = VoxelizedModelImageType::New();
  label->CopyInformation( inputImage );
  label->SetRegions( label->GetLargestPossibleRegion() );
  label->Allocate();
  label->FillBuffer( 0 );

  //Calculation the sample-spacing, i.e the half of the smallest spacing existing in the original image
  double minSpacing = label->GetSpacing()[0];
  for (unsigned int i = 1; i < label->GetSpacing().Size(); ++i)
    {
    minSpacing = std::min(label->GetSpacing()[i], minSpacing);
    }

  vtkNew<vtkPolyDataPointIDSampler> sampler;

  sampler->SetInput( model );
  sampler->SetDistance( samplingDistance * minSpacing );
  sampler->GenerateEdgePointsOn();
  sampler->GenerateInteriorPointsOn();
  sampler->GenerateVertexPointsOn();
  sampler->Update();

  std::cout << model->GetNumberOfPoints() << std::endl;
  std::cout << sampler->GetOutput()->GetNumberOfPoints() << std::endl;

  vtkPoints* points = sampler->GetOutput()->GetPoints();
  vtkPointData* pointData = sampler->GetOutput()->GetPointData();
  vtkIdTypeArray* modelPointIndexes = vtkIdTypeArray::SafeDownCast(
    pointData->GetScalars("pointIndexes"));
  for( vtkIdType k = 0; k < points->GetNumberOfPoints(); k++ )
    {
    double* pt = points->GetPoint( k );
    VoxelizedModelImageType::PointType pitk;
    pitk[0] = pt[0];
    pitk[1] = pt[1];
    pitk[2] = pt[2];
    VoxelizedModelImageType::IndexType idx;
    label->TransformPhysicalPointToIndex( pitk, idx );

    if( label->GetLargestPossibleRegion().IsInside(idx) )
      {
      vtkIdType pointId = modelPointIndexes->GetValue(k);
      label->SetPixel( idx, pointId );
      }
    }
  return label;
}

//-----------------------------------------------------------------------------
vtkPolyData* ReadPolyData(const std::string& fileName)
{
  vtkPolyData* polyData = 0;
  vtkSmartPointer<vtkPolyDataReader> pdReader;
  vtkSmartPointer<vtkXMLPolyDataReader> pdxReader;

  // do we have vtk or vtp models?
  std::string::size_type loc = fileName.find_last_of(".");
  if( loc == std::string::npos )
    {
    std::cerr << "Failed to find an extension for " << fileName << std::endl;
    return polyData;
    }

  std::string extension = fileName.substr(loc);

  if( extension == std::string(".vtk") )
    {
    pdReader = vtkSmartPointer<vtkPolyDataReader>::New();
    pdReader->SetFileName(fileName.c_str() );
    pdReader->Update();
    polyData = pdReader->GetOutput();
    }
  else if( extension == std::string(".vtp") )
    {
    pdxReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    pdxReader->SetFileName(fileName.c_str() );
    pdxReader->Update();
    polyData = pdxReader->GetOutput();
    }
  if( polyData == NULL )
    {
    std::cerr << "Failed to read surface " << fileName << std::endl;
    return polyData;
    }

  // LPS vs RAS

  vtkPoints * allPoints = polyData->GetPoints();
  for( int k = 0; k < allPoints->GetNumberOfPoints(); k++ )
    {
    double* point = polyData->GetPoint( k );
    point[0] = -point[0];
    point[1] = -point[1];
    allPoints->SetPoint( k, point[0], point[1], point[2] );
    }
  polyData->Register(0);
  return polyData;
}

//-----------------------------------------------------------------------------
template <class T>
void WriteImage(typename InputImage<T>::Type::Pointer image,
                const std::string& fileName,
                ModuleProcessInformation* processInformation)
{
  typedef itk::ImageFileWriter<InputImage<T>::Type> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  itk::PluginFilterWatcher watchWriter(writer,
                                       "Write Volume",
                                       processInformation);
  writer->SetFileName( fileName.c_str() );
  writer->SetInput( image );
  writer->SetUseCompression(1);
  writer->Update();
}
