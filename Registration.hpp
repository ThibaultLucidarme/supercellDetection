
#ifndef _Registration_hpp
#define _Registration_hpp

#include <iostream>
#include "Types.hpp"

#include <itkImageRegistrationMethod.h>

#include <itkTranslationTransform.h>
#include <itkEuler3DTransform.h>
#include <itkRigid3DTransform.h>

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMutualInformationImageToImageMetric.h>

#include <itkLinearInterpolateImageFunction.h>

#include <itkRegularStepGradientDescentOptimizer.h>

#include <itkResampleImageFilter.h>




#include <itkCommand.h>

namespace itk
{
    
    class CommandIteration : public itk::Command {
    public:
        typedef CommandIteration  Self;
        typedef itk::Command         SuperClass;
        typedef itk::SmartPointer< Self >    Pointer;
      itkNewMacro( Self );
    protected:
        CommandIteration()
        {
        }
        
    public:
        
        typedef itk::RegularStepGradientDescentOptimizer   OptimizerType;
        typedef const OptimizerType  * OptimizerPointer;
        
        void Execute(  itk::Object * caller, const itk::EventObject & event )
        {
            this-> Execute(  (const itk::Object *) caller,  event  );
        }
        
        void Execute(  const itk::Object * caller, const itk::EventObject & event )
      {
	OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( caller );
            if( typeid( event  ) == typeid( itk::IterationEvent ) )
            {
                std::cout << optimizer->GetCurrentIteration() << " : ";
                std::cout << optimizer->GetValue() << " : ";
                std::cout << optimizer->GetCurrentPosition() << std::endl;
            }
        }
    };
    
    
    template< class T>
    class RegistrationFilter : public ImageToImageFilter< T, T >
    {
        
    private:
        
        typedef itk::Euler3DTransform< double > TransformType;
        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        typedef itk::LinearInterpolateImageFunction< ImageType , double > InterpolatorType;
        //typedef itk::MeanSquaresImageToImageMetric< ImageType , ImageType > MetricType;
        typedef itk::MutualInformationImageToImageMetric< ImageType , ImageType > MetricType;
        typedef itk::ImageRegistrationMethod< ImageType , ImageType > RegistrationType;
        typedef RegistrationType::ParametersType     ParametersType;
        typedef itk::ResampleImageFilter< ImageType , ImageType > ResamplerType;
        
        
        
    public:
        /** Standard class typedefs. */
        typedef RegistrationFilter			 Self;
        typedef ImageToImageFilter< T, T >	 Superclass;
        typedef SmartPointer< Self >		 Pointer;
        typedef SmartPointer< const Self >	 ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(RegistrationFilter, ImageToImageFilter);
        
        /** Parameter personalization */
        void SetFixedImage(typename T::Pointer fi)
        {
       	  this->ProcessObject::SetNthInput( 0,  fi );
        }
        void SetMovingImage(typename T::Pointer mi)
        {
	  this->ProcessObject::SetNthInput( 1, mi );
        }
        void ObserveRegistration(bool obr)
        {
            _observeRegistration = obr;
        }
        
    protected:
        RegistrationFilter()
        {
	  this->SetNumberOfRequiredInputs(2);

	  
            _transform      = TransformType::New();
            _transform->SetIdentity();
            _optimizer      = OptimizerType::New();
            _optimizer->SetMaximumStepLength( 4.00 );
            _optimizer->SetMinimumStepLength( 0.01  );
            _optimizer->SetNumberOfIterations( 100  );
            _optimizer->MaximizeOff();
            
            _observeRegistration = true;
            _observer  = CommandIteration::New();
            
            _interpolator   = InterpolatorType::New();
            
            _metric         = MetricType::New();
            
            
            _registrator    = RegistrationType::New();
            _registrator->SetTransform(     _transform    );
            _registrator->SetOptimizer(     _optimizer    );
            _registrator->SetInterpolator(  _interpolator );
            _registrator->SetMetric(        _metric       );
            
            _registrator->SetInitialTransformParameters( _transform->GetParameters() );
            
            _resampler  = ResamplerType::New();
            _resampler->SetTransform ( _transform );
            
            
        }
        
        // Does the real work.
        virtual void GenerateData()
        {
            std::cout<<" --- BEGIN REGISTRATION --- "<<std::endl;
            
            std::cout<<"Setting parameters"<<std::endl;
       	  _fixedImage = const_cast<T*>( this->GetInput(0) );
          _movingImage = const_cast<T*>( this->GetInput(1) );
            _registrator->SetFixedImage(    _fixedImage );
            _registrator->SetMovingImage(   _movingImage );

	    _fixedImage->Print(std::cout);

            _registrator->SetFixedImageRegion( _fixedImage->GetLargestPossibleRegion() );

            _resampler->SetInput( _movingImage );
            _resampler->SetOutputOrigin( _fixedImage->GetOrigin() );
            _resampler->SetOutputSpacing( _fixedImage->GetSpacing() );
            _resampler->SetSize( _fixedImage->GetLargestPossibleRegion().GetSize() );
            
            if ( _observeRegistration )
                _optimizer->AddObserver( itk::IterationEvent(), _observer);
            
            try
            {
	      std::cout<<"Optimizing transform"<<std::endl;
                _registrator->Update();
            }
            catch( itk::ExceptionObject & excp )
            {
                std::cerr << "Error in registration" << std::endl;
                std::cerr << excp << std::endl;
            }
            
            _transform->SetParameters( _registrator->GetLastTransformParameters()  );
            
            std::cout<<"Resampling moving image from optimized transform"<<std::endl;
            //_resampler->Update();
            
            std::cout<<" --- END REGISTRATION --- "<<std::endl;
            
            this->GraftOutput( _fixedImage );
            
        }
        
        RegistrationFilter(const Self &);  //purposely not implemented
        void operator=(const Self &);    //purposely not implemented
        
        //Component filters
        typename TransformType::Pointer     _transform;
        typename OptimizerType::Pointer     _optimizer;
        typename InterpolatorType::Pointer  _interpolator;
        typename MetricType::Pointer        _metric;
        typename RegistrationType::Pointer  _registrator;
        typename ResamplerType::Pointer     _resampler;
        typename CommandIteration::Pointer  _observer;
        
        typename T::Pointer     _fixedImage, _movingImage;
        
        bool _observeRegistration;
        
    };
    
} //namespace ITK



#endif /* defined(____Registration__) */
