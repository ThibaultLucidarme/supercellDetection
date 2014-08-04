#ifndef _main_hpp
#define _main_hpp

//ITK includes
#include <itkImage.h>
#include <itkExceptionObject.h>
#include <itkSimpleFilterWatcher.h>
#include <itkRGBPixel.h>
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"
#include "itkRandomImageSource.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkImageDuplicator.h"
//STD includes
#include <iostream>
#include <string>

//Project includes
#include "Types.hpp"
#include "CommandLineParser.hpp"

/*
 ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
 progress->SetMiniPipelineFilter(this);
 progress->RegisterInternalFilter(erode, .5f);
 progress->RegisterInternalFilter(dilate, .5f);
 */

#endif
