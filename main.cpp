#include "main.hpp"
#include "Registration.hpp"


int main(int argc, char *argv[])
{

    /***************************************
     UTILITY
	Self implemented alternative to the boost command line parser
     **************************************/
	 
    p::CommandLineParser parser(argc, argv);
    
    bool verbose = parser.addOption<bool>("-v", false, "display all progress details");
    std::string dirName = parser.addOption<std::string>("-d","/home/tlucidarme/Workspace/dataset2_proj/","directory contaning the input files (has to end with '/')");
    unsigned int imageNum = parser.addOption<unsigned int>("-n",30,"image number to read in directory provided with -d option");
    
    typedef itk::ImageFileReader< RGBImageType >	ReaderType;
    typedef itk::RGBToLuminanceImageFilter< RGBImageType, ImageType > RGB2GRAYFilterType;
    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    
    /***************************************
     Read data dir
     **************************************/
    
    //char tmp[7];
    //sprintf (tmp,"proj_%05d.bmp", imageNum);
    //std::string movingimagefilename = dirName+tmp;
    std::string movingimagefilename = dirName+"proj_raw_00030.bmp";
    
	//read
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(movingimagefilename);
    
	//convert to gray
    RGB2GRAYFilterType::Pointer toGrayFilter = RGB2GRAYFilterType::New();
    toGrayFilter->SetInput( reader->GetOutput() );
    toGrayFilter->Update();// needed for saving movingImage
	
	//save resulting image
    ImageType::Pointer movingImage = toGrayFilter->GetOutput();
    
    /***************************************
     Read each pair of images (moving and fixed)
     **************************************/
    
    
    for(int frame = 180 ; frame<10300; frame+=150 )
    {
        // previous moving image becomes fixed
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(movingImage);
        duplicator->Update();
        ImageType::Pointer fixedImage = duplicator->GetOutput();
        std::string fixedimagefilename = movingimagefilename;
        
        // Get new moving image
        sprintf (tmp,"proj_raw_%05d.bmp", frame);
        movingimagefilename = dirName+tmp;
        
        reader = ReaderType::New();
        reader->SetFileName(movingimagefilename);
        reader->UpdateLargestPossibleRegion();
        
        toGrayFilter = RGB2GRAYFilterType::New();
        toGrayFilter->SetInput( reader->GetOutput() );
        toGrayFilter->Update();
        
        movingImage = toGrayFilter->GetOutput();
        
        std::cerr<< "Registering "<< fixedimagefilename <<" and "<< movingimagefilename <<std::endl;
        
        /***************************************
         Registration (rewritten overall in RegistrationFilter class)
         **************************************/
        
		//-------
		// Type definitions:
		// controls which algorithms are used within the registration process
		//-------
		
        //typedef itk::Euler3DTransform< double > TransformType;
        //typedef itk::TranslationTransform< double, Dimension> TransformType;
        typedef itk::Rigid3DTransform< double> TransformType;
        
        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        
        typedef itk::LinearInterpolateImageFunction< ImageType , double > InterpolatorType;
        
        typedef itk::MutualInformationImageToImageMetric< ImageType , ImageType > MetricType;
        //typedef itk::MeanSquaresImageToImageMetric< ImageType , ImageType > MetricType;
        
        typedef itk::MultiResolutionImageRegistrationMethod< ImageType , ImageType > RegistrationType;
        
        typedef RegistrationType::ParametersType     ParametersType;
        
        typedef itk::ResampleImageFilter< ImageType , ImageType > ResamplerType;
        
        typedef itk::MultiResolutionPyramidImageFilter< ImageType , ImageType > FixedImagePyramidType;
        typedef itk::MultiResolutionPyramidImageFilter< ImageType , ImageType > MovingImagePyramidType;
        
		//-------
		// Initialization
		//-------
        TransformType::Pointer      transform     = TransformType::New();
        OptimizerType::Pointer      optimizer     = OptimizerType::New();
        InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
        RegistrationType::Pointer   registration  = RegistrationType::New();
        MetricType::Pointer         metric        = MetricType::New();
        FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
        MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();
        
		//-------
		// Instanciation and pipelining
		//-------
        transform = TransformType::New();
        transform->SetIdentity();
        
        optimizer = OptimizerType::New();
        optimizer->SetMaximumStepLength( 4.00 );
        optimizer->SetMinimumStepLength( 0.01  );
        optimizer->SetNumberOfIterations( 100  );
        optimizer->MaximizeOff();
        optimizer->SetGradientMagnitudeTolerance(0.00001);
        
        itk::CommandIteration::Pointer observer;
        observer  = itk::CommandIteration::New();
        if ( verbose )
            optimizer->AddObserver( itk::IterationEvent(), observer);
        
        
        interpolator   = InterpolatorType::New();
        
        metric         = MetricType::New();
        
        registration->SetOptimizer(     optimizer     );
        //registration->SetInitialTransformParameters( transform->GetParameters() );
        registration->SetTransform(     transform     );
        registration->SetInterpolator(  interpolator  );
        registration->SetMetric( metric  );
        registration->SetFixedImagePyramid( fixedImagePyramid );
        registration->SetMovingImagePyramid( movingImagePyramid );
        
        registration->SetFixedImage(  fixedImage  );
        registration->SetMovingImage(   movingImage   );
        registration->SetFixedImageRegion( fixedImage->GetBufferedRegion() );
        registration->SetNumberOfLevels( 3 );
        
        ParametersType initialParameters( transform->GetNumberOfParameters() );
        for (int i=0; i<transform->GetNumberOfParameters(); i++)
            initialParameters[i] = 0.0;
        registration->SetInitialTransformParameters( initialParameters );
        
        
        // */
        
        /***************************************
         Launch Registration Pipeline to get transform
         ***************************************/
        if(verbose)
            itk::SimpleFilterWatcher watcher(registration, "registration");
        
        try
        {
            registration->Update();
            std::cout << ""
            << registration->GetOptimizer()->GetStopConditionDescription()
            << std::endl;
        }
        catch( itk::ExceptionObject & err )
        {
            std::cout << "ExceptionObject caught !" << std::endl;
            std::cout << err << std::endl;
            return EXIT_FAILURE;
        }
		
		// */
		
		/***************************************
         Apply transform to moving image through resampler
         ***************************************/
        
        typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;
        TransformType::Pointer finalTransform = TransformType::New();
        ParametersType finalParameters = registration->GetLastTransformParameters();
        
        finalTransform->SetParameters( finalParameters );
        finalTransform->SetFixedParameters( transform->GetFixedParameters() );
        
        ResampleFilterType::Pointer resample = ResampleFilterType::New();
        
        resample->SetTransform( finalTransform );
        resample->SetInput( movingImage );
        
        resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
        resample->SetOutputOrigin(  fixedImage->GetOrigin() );
        resample->SetOutputSpacing( fixedImage->GetSpacing() );
        resample->SetOutputDirection( fixedImage->GetDirection() );
        resample->SetDefaultPixelValue( 100.0 );
        
        resample->Update();
        
		// */
		
        /***************************************
         Write resulting Images
         **************************************/
         
        typedef itk::ImageFileWriter< ImageType > WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput(resample->GetOutput());
        sprintf (tmp,"proj_raw_%05d_registered.bmp", frame);
        std::cout<<"##########"<<dirName+tmp<<std::endl;
        writer->SetFileName(dirName+tmp);
        writer->UpdateLargestPossibleRegion();
        
         
        // */
    }
    
    std::cout<<"\nSUCESS"<<std::endl;
    
    return EXIT_SUCCESS;
    
}
