
#ifndef _itkWatershed_hpp
#define _itkWatershed_hpp

#include <iostream>
#include "Types.hpp"

#include <itkWatershedImageFilter.h>
#include <itkCastImageFilter.h>

namespace itk
{
    template< class T>
    class WatershedFilter : public ImageToImageFilter< T, T >
    {
        
    private:
        
        typedef itk::WatershedImageFilter<T> WatershedFilterType;
        typedef itk::Image<unsigned long, Dimension> ImageLongType;
        typedef itk::CastImageFilter< ImageLongType, T > CastFilterType;
        
        
    public:
        /** Standard class typedefs. */
        typedef WatershedFilter				 Self;
        typedef ImageToImageFilter< T, T >	 Superclass;
        typedef SmartPointer< Self >		 Pointer;
        typedef SmartPointer< const Self >	 ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(WatershedFilter, ImageToImageFilter);
        
        /** Parameter personalization */
        void SetLevel(float lvl)
        {
            _level = lvl;
        }
        void SetLowerThesh(float lt)
        {
            _lthresh =  lt;
        }
        
    protected:
        WatershedFilter()
        {
            //typename T::ConstPointer input = this->GetInput();
            //typename T::Pointer output = this->GetOutput();
            
            /*
             Finally we set up the watershed filter. There are two parameters. Level controls watershed depth, and Threshold controls the lower thresholding of the input. Both parameters are set as a percentage (0.0 - 1.0) of the maximum depth in the input image.
             */
            _watershed = WatershedFilterType::New();
            _watershed->SetLevel(_level);
            _watershed->SetThreshold(_lthresh);
            
            _caster = CastFilterType::New();
            
        }
        
        // Does the real work.
        virtual void GenerateData()
        {
            _watershed->SetInput( this->GetInput() );
            _caster->SetInput( _watershed->GetOutput() );
            
            _caster->Update();
            
            this->GraftOutput( _caster->GetOutput() );
        }
     
        WatershedFilter(const Self &);  //purposely not implemented
        void operator=(const Self &);    //purposely not implemented
        void PrintSelf(std::ostream& os, itk::Indent indent) const
        {
            Superclass::PrintSelf(os, indent);
            
            os << indent << "level:" << this->_level << std::endl;
            os << indent << "Lower threshold:" << this->_lthresh << std::endl;
        }
        
        //Component filters
        typename CastFilterType::Pointer _caster;
        typename WatershedFilterType::Pointer _watershed;
        
        
        float _level;
        float _lthresh;
    
    };
    
} //namespace ITK



#endif /* defined(____Watershed__) */
