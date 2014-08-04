
#ifndef _Types_hpp
#define _Types_hpp

#include <itkImage.h>


const unsigned int Dimension = 2;

/*************global typedefs**********/

//Pixel types
typedef double                                  PixelType;
typedef itk::RGBPixel<unsigned char>            RGBPixelType;

//Image types
typedef itk::Image< PixelType, Dimension>		ImageType;
typedef itk::Image< unsigned long, Dimension>   ImageLabelType;
typedef itk::Image< RGBPixelType, Dimension>	RGBImageType;

#endif
