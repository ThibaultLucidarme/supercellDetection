#ifndef _Reader_hpp
#define _Reader_hpp

//#include <netcdf.h>
#include <iostream>
#include <exception>

#include "/home/tlucidarme/Library/Developer/NetCDFcpp/build/include/netcdfcpp.h"

#include "Types.hpp"

#include "itkImageToImageFilter.h"
#include "itkImportImageFilter.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}

    /***************************************
     Reads NetCDF files and return ITK images
	 TODO:
	 transform into Factory to support HDF4, HDF5, DICOM image readers
     **************************************/

template< class T>
class Reader
{
public:
    
    // Init importFilter
    Reader()
    {
        _importFilter = ImportFilterType::New();
    }
    
    ~Reader()
    {
        delete[] _data;
    }
    
    // Create itk Image::Pointer
    void Update()
    {
        //get data
        ReadFile();
        
        //process
        _importFilter->Update();
    }
    
    /** Parameter personalization */
    void SetInputFile(std::string f)
    {
        _filename = f;
    }
    void SetVariable(std::string v)
    {
        _varname = v;
    }
    ImageType::Pointer GetOutput()
    {
        return _importFilter->GetOutput();
    }
    
private:
    
    
    // Read netCDF
    void Read()
    {
        
        // Open the file.
        NcFile file( _filename.c_str(), NcFile::ReadOnly);
        
        NcVar* var = file.get_var( _varname.c_str() );
        
        
        if( var->is_valid() )
        {
            _sizeFlat = var->num_vals();
            _sizeArray = var->edges();
            
            /*
             _dim = var->num_dims();
             
             std::cout<< _dim<<std::endl;
             for (int i=0; i<_dim; i++)
             std::cout<<_sizeArray[i]<<"\t";
             std::cout<<"\t"<<_sizeFlat<<std::endl;
             // */
            
            _data = new PixelType[_sizeFlat];
            
            try
            {
                var->get(_data, _sizeArray);
            }
            catch(std::exception& a)
            {
                a.what();
                
                //free mem
                delete[] _data;
                _sizeArray = NULL;
                var = NULL;
                
                //quit
                std::cout<<"Could not read variable :"<<_varname<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            
        }
        else
        {
            var = NULL;
            
            //quit
            std::cout<<"Could not find variable :"<<_varname<<std::endl;
            exit(EXIT_FAILURE);
        }
        
    }
    
    // Create itk Image::Pointer from buffer
    void ReadFile()
    {
        //fill _data
        Read();
        
        ImportFilterType::SizeType  size;
        
        size[0]  = _sizeArray[1];  // size along X
        size[1]  = _sizeArray[2];  // size along Y
        size[2]  = _sizeArray[3];  // size along Z
        
        ImportFilterType::IndexType start;
        start.Fill( 0 );
        
        ImportFilterType::RegionType region;
        region.SetIndex( start );
        region.SetSize(  size  );
        
        _importFilter->SetRegion( region );
        
        double origin[ Dimension ];
        origin[0] = 0.0;    // X coordinate
        origin[1] = 0.0;    // Y coordinate
        origin[2] = 0.0;    // Z coordinate
        
        _importFilter->SetOrigin( origin );
        
        double spacing[ Dimension ];
        spacing[0] = 1.0;    // along X direction
        spacing[1] = 1.0;    // along Y direction
        spacing[2] = 1.0;    // along Z direction
        
        _importFilter->SetSpacing( spacing );
        
        const bool importImageFilterWillOwnTheBuffer = true;
        _importFilter->SetImportPointer( _data, _sizeFlat,
                                        importImageFilterWillOwnTheBuffer );
    }
    
protected:
    
    typedef itk::ImportImageFilter< PixelType, Dimension >   ImportFilterType;
    
    //Component filters
    ImportFilterType::Pointer _importFilter;
    
    std::string _filename, _varname;
    int _dim;
    long* _sizeArray;
    long _sizeFlat;
    PixelType* _data;
    ImageType::Pointer _img;
};


#endif
