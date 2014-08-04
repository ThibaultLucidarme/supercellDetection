#ifndef _Writer_hpp
#define _Writer_hpp

//#include <netcdf.h>
#include <iostream>
#include <exception>

#include "/home/tlucidarme/Library/Developer/NetCDFcpp/build/include/netcdfcpp.h"

#include "Types.hpp"
#include <itkImage.h>


    /***************************************
     Writes NetCDF files and return ITK images
	 TODO:
	 transform into Factory to support HDF4, HDF5, DICOM image writers
     **************************************/

template< class T>
class Writer
{
public:
    
    // Init importFilter
    Writer()
    {
    }
    
    ~Writer()
    {
    }
    
    // Create itk Image::Pointer
    void Update()
    {
      std::cout<<"Writing "<<_filename<<std::endl;
        //get data
        WriteFile();
    }
    
    /** Parameter personalization */
    void SetOutputFile(std::string f)
    {
        _filename = f;
    }
    void SetVariable(std::string v)
    {
        _varname = v;
    }
    void SetInput(typename T::Pointer i)
    {
        _img = i;
    }
    
private:
    
    
    // Read netCDF
    void WriteFile()
    {   
        //Open the file.
        NcFile file( _filename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
        
        //set dimension
        typename T::SizeType size = _img->GetLargestPossibleRegion().GetSize();
        
        // Define the dimensions. NetCDF will hand back an ncDim object for each.
        NcDim * lvlDim = file.add_dim("level", size[0] );
        NcDim * latDim = file.add_dim("latitude", size[1] );
        NcDim * lonDim = file.add_dim("longitude", size[2] );
        NcDim * recDim = file.add_dim("record", 1 );  //adds an unlimited dimension
        
        //Define the netCDF variable
        NcVar * var = file.add_var( _varname.c_str(), ncFloat, recDim, lvlDim, latDim, lonDim);
        
        //write data
        var->put_rec( _img->GetBufferPointer(), 0 );
        
        //free mem
        file.close();
        
    }
    
    
protected:
    
    std::string _filename, _varname;
    typename T::Pointer _img;
};


#endif
