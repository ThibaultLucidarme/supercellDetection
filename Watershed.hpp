
#ifndef _lakWatershed_hpp
#define _lakWatershed_hpp

    /***************************************
     Self implementation of 3D Watershed algorithm to replace default ITK one
	 implemented following the following article:
	 
@Article{stormattr,
author = {Valliappa Lakshmanan and Travis Smith},
title = {Data Mining Storm Attributes from Spatial Grids},
journal = {J. Ocea. and Atmos. Tech.},
year = {2009},
volume = {26},
number = {11},
pages = {2353-2365}
}

3D data is stored as flattened 1D array

     **************************************/

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <set>
#include <float.h>
//#include <boost/unordered_set.hpp>

#include "Types.hpp"

#include "itkImageToImageFilter.h"
#include <itkImageIterator.h>

#include "itkWatershedImageFilter.h"
#include <itkCastImageFilter.h>

    /***************************************
     Levelset: Manages the various levelsets for the flooding (topographic representation)
     **************************************/

class LevelSet
{
public:

    // Input
    std::vector<PixelType> gridVals;
    
    // Output
    std::vector< std::vector<int> > binIdcs;
    std::vector<PixelType> binVals;
    
    // Params
    double stepSize;
    PixelType lThresh;
    PixelType uThresh;
    
    int nBins;
    
    // Clock for timing
    std::clock_t start;
    
    // METHODS
    
    /*
     // Construtors / destructors
     LevelSet(int);
     LevelSet(std::vector<double> *, double, double, double);
     // Don't think I need to do anything w/ dynamic memory right now, so no destructor yet.
     
     
     // Methods
     void runAll();
     void memoryEfficient();
     void equalDistAssumption();
     void dontRecalcBinIdx();
     */
    
    LevelSet(int size)
    {
        // Set parameters
        stepSize = 2.0;
        lThresh = 20.0; // Inclusive threshold
        uThresh = 70.0; // Exclusive threshold
        
        // Create "random" variable
        gridVals.resize(size);
        for (size_t idx=0; idx < gridVals.size(); idx++)
        {
            //        gridVals[idx] = 100.0 * (double) (rand() + 1) / (double) (RAND_MAX + 1);
            gridVals[idx] = 20.0 + 50.0 * (PixelType) (rand() + 1) / (PixelType) (RAND_MAX + 1);
        }
        
    }
    
    LevelSet(std::vector<PixelType> * lGridVals, double lStepSize, PixelType lLThresh, PixelType lUThresh)
    {
        gridVals = * lGridVals;
        stepSize = lStepSize;
        lThresh = lLThresh;
        uThresh = lUThresh;
        
        memoryEfficient();
    }
    
    void runAll()
    {
        // Just so resize doesn't mess w/ timing...
        dontRecalcBinIdx();
        
        std::cout << "\n\nBEGIN MAIN RUNS\n\n" << std::endl;
        start = std::clock();
        memoryEfficient();
        std::cout << (std::clock() - start) / ((double) CLOCKS_PER_SEC) << std::endl;
        
        start = std::clock();
        equalDistAssumption();
        std::cout << (std::clock() - start) / ((double) CLOCKS_PER_SEC) << std::endl;
        
        start = std::clock();
        dontRecalcBinIdx();
        std::cout << (std::clock() - start) / ((double) CLOCKS_PER_SEC) << std::endl;
        
        return;
    }
    
	// Improve memory efficiency using preallocation
    void memoryEfficient()
    {
        //std::cout << "BEGIN memoryEfficient" << std::endl;
        
        //std::cout << "A" << std::endl;
        // Determine number of elements in each bin
        nBins = (int) ceil((uThresh - lThresh) / stepSize);
        std::vector<int> binSizes(nBins, 0);
        for (size_t idx=0; idx < gridVals.size(); idx++)
        {
            if ( (gridVals[idx] >= lThresh) && (gridVals[idx] < uThresh) )
            {
                int binIdx = (int) floor((gridVals[idx] - lThresh) / stepSize);
                binSizes[binIdx] += 1;
            }
        }
        
        //std::cout << "B" << std::endl;
        // Allocate appropriate amount of memory in binIdcs and binVals
        binIdcs.resize(nBins);
        binVals.resize(nBins);
        for (int binIdx=0; binIdx < nBins; binIdx++)
        {
            double binVal = lThresh + ((double) binIdx) * stepSize;
            binVals[binIdx] = binVal;
            
            std::vector<int> idcs;
            idcs.reserve(binSizes[binIdx]);
            binIdcs[binIdx] = idcs;
        }
        
        //std::cout << "C" << std::endl;
        // Assign indicies to appropirate bin
        for (size_t idx=0; idx < gridVals.size(); idx++)
        {
            if ( (gridVals[idx] >= lThresh) && (gridVals[idx] < uThresh) )
            {
                int binIdx = (int) floor((gridVals[idx] - lThresh) / stepSize);
                binIdcs[binIdx].push_back(idx);
            }
        }
        
        return;
    }
    
    //Assumes relative equalisty in gray level distribution
    void equalDistAssumption()
    {
        std::cout << "BEGIN equalDistAssumption" << std::endl;
        
        
        // This test somewhat broken by using the random values ...
        
        std::cout << "A" << std::endl;
        // Allocate appropriate amount of memory in binIdcs, binVals, and lastIdx
        nBins = (int) ceil((uThresh - lThresh) / stepSize);
        binIdcs.resize(nBins);
        binVals.resize(nBins);
        std::vector<int> lastIdx(nBins, 0);
        for (int binIdx=0; binIdx < nBins; binIdx++)
        {
            PixelType binVal = lThresh + ((PixelType) binIdx) * stepSize;
            binVals[binIdx] = binVal;
            
            std::vector<int> idcs;
            idcs.reserve(gridVals.size() / nBins); // Assume initial distribution of values is equally divided among the bins
            binIdcs[binIdx] = idcs;
        }
        
        std::cout << "B" << std::endl;
        // Assign grid indicies to appropriate bin
        for (size_t idx=0; idx < gridVals.size(); idx++)
        {
            int binIdx = (int) floor((gridVals[idx] - lThresh) / stepSize);
            if ( (binIdx >= 0) && (binIdx < nBins) )
            {
                binIdcs[binIdx].push_back(idx);
                lastIdx[binIdx] += 1;
                if (lastIdx[binIdx] >= (int) binIdcs[binIdx].size())
                {
                    binIdcs[binIdx].resize((int) ((PixelType) binIdcs[binIdx].size() + ((PixelType) gridVals.size() / (PixelType) nBins)));
                }
            }
        }
        
        return;
    }
    
    void dontRecalcBinIdx()
    {
        std::cout << "BEGIN dontRecalcBinIdx" << std::endl;
        
        std::cout << "A" << std::endl;
        // Determine number of elements in each bin
        nBins = (int) ceil((uThresh - lThresh) / stepSize);
        std::vector<int> gridBinIdx(gridVals.size(), nBins);
        std::vector<int> binSizes(nBins, 0);
        for (size_t idx=0; idx < gridVals.size(); idx++)
        {
            if ( (gridVals[idx] >= lThresh) && (gridVals[idx] < uThresh) )
            {
                int binIdx = (int) floor((gridVals[idx] - lThresh) / stepSize);
                binSizes[binIdx] += 1;
                gridBinIdx[idx] = binIdx;
            }
        }
        
        std::cout << "B" << std::endl;
        // Allocate appropriate amount of memory in binIdcs, binVals, and lastIdx
        binIdcs.resize(nBins);
        binVals.resize(nBins);
        std::vector<int> lastIdx(nBins, 0);
        for (int binIdx=0; binIdx < nBins; binIdx++)
        {
            PixelType binVal = lThresh + ((PixelType) binIdx) * stepSize;
            binVals[binIdx] = binVal;
            
            std::vector<int> idcs(binSizes[binIdx]);
            binIdcs[binIdx] = idcs;
        }
        
        std::cout << "C" << std::endl;
        // Assign indicies to appropriate bin
        for (size_t idx=0; idx < gridBinIdx.size(); idx++)
        {
            if (gridBinIdx[idx] != nBins)
            {
                int binIdx = gridBinIdx[idx];
                binIdcs[binIdx][lastIdx[binIdx]] = idx;
                lastIdx[binIdx] += 1;
            }
        }
        
        return;
    }
    
};

    /***************************************
     Label: Manages the labels assigned to each flood basins
     **************************************/

class Label
{
public:
    // ATTRIBUTES
    int value;      // Input value of this label
    int pixelCount; // Number of pixels with this label
    
    const static int MIN_PIXEL_COUNT = 300; // Mininum number of pixels a label is required to have before it cannot be merged with larger labels
    
    static int numLabels;   // Counter of number of labels in existence. Increases on Label creation, decreases on Label merger, NOT on destruction
    
    // METHODS
    /*
     Label(int); // Constructor
     
     void mergeWithLabel(Label *);   // Merges this label with input label
     void addToPixelCount(int);      // Increases pixelCount by input value
     
     friend bool operator== (Label &l1, Label &l2);  // == operator override
     friend bool operator!= (Label &l1, Label &l2);  // != operator override
     */
    
    Label(int lValue)
    {
        pixelCount = 0; // Set pixelCount
        value = lValue; // Set value
        
        numLabels++;    // Increase static number of labels
    }
    
    void addToPixelCount(int num)
    {
        pixelCount += num;
    }
    
    
    void mergeWithLabel(Label * mLabel)
    {
        mLabel->addToPixelCount(pixelCount);    // Add this label's pixelCount to other label's
        pixelCount = mLabel->pixelCount;        // Set this label's pixelCount to other label's
        value = mLabel->value;                  // Set this label's value to other label's
        
        numLabels--;    // Decrease static number of labels
    }
    
    
    
};

bool operator==(Label &l1, Label &l2)
{
    return ( (l1.value == l2.value) && (l1.pixelCount == l2.pixelCount) );
}

bool operator!=(Label &l1, Label &l2)
{
    return !(l1 == l2);
}

    /***************************************
     Pixel: manages label and value assigned to the image during the region growth
     **************************************/

class Pixel
{
public:
    // ATTRIBUTES
    int idx;        // Input index of pixel on grid
    bool onEdge;    // Is the pixel on the grid edge
    Label * label;  // Pointer to label of the pixel
    
    std::vector<Label *> neighborLabels;    // Pointers to labels of neighboring pixels
    
    static std::vector<int> neighborMod;
    
    const static int NUM_NEIGHBORS = 26;
    
    // METHODS
    /*
     Pixel(int,int,int,int,Label *); // Constructor
     
     void addNeighborLabel(Label *); // Add pointer to label of neighboring pixel
     bool tryToLabel();              // Try to set an appropriate label value for this pixel - return true if is set, false if not set or already labeled
     bool isLabeled();               // True if pixel has been labeled with a valid value
     
     static void setupNeighborMod(int,int,int);  // Setup neighborMod for current run
     */
    
    
    Pixel(int pIdx, int iSize, int jSize, int kSize, Label * pLabel)
    {
        neighborLabels.reserve(NUM_NEIGHBORS);   // Allocate space for pointers to neighboring labels
        
        label = pLabel; // Set label pointer
        
        idx = pIdx; // Set index value;
        
        // Get i, j, k values
        int tmpIdx = idx+1;
        int k = tmpIdx / (iSize * jSize);
        
        tmpIdx -= k * iSize * jSize;
        int j = tmpIdx / iSize;
        
        tmpIdx -= j * iSize;
        int i = tmpIdx - 1;
        
        // Check if index is on edge of grid
        if ( (i == 0) || (j == 0) || (k == 0) || (i == iSize - 1) || (j == jSize - 1) || (k == kSize - 1) )
        {
            onEdge = true;
        }
        else
        {
            onEdge = false;
        }
    }
    
    
	// get 3x3x3 neighborhood of the 1D image
    static void setupNeighborMod(int iSize, int jSize, int kSize)
    {
        neighborMod.reserve(Pixel::NUM_NEIGHBORS);
        for (int i=-1; i <= 1; i++)
            for (int j=-1; j <= 1; j++)
                for (int k=-1; k <= 1; k++)
                    if ( (i != 0) && (j != 0) && (k != 0) )
                    {
                        int modVal = i + iSize * j + iSize * jSize * k -1;
                        neighborMod.push_back(modVal);
                    }

    }
    
    void addNeighborLabel(Label * neighborLabel)
    {
        neighborLabels.push_back(neighborLabel);
    }
    
    bool tryToLabel()
    {
        // Pixel is on edge and not considered
        if (onEdge)
        {
            //        std::cout << "tryToLabel: A" << std::endl;
            return false;
        }
        // Pixel is already labeled
        else if (isLabeled())
        {
            //        std::cout << "tryToLabel: B" << std::endl;
            return false;
        }
        // Pixel needs to be labeled
        else
        {
            //        std::cout << "tryToLabel: C" << std::endl;
            // Using vector - unordered_set got too wonky; access to correct memory region uncertain
            // Loop through neighbor labels
            std::vector<Label *> uniqueLabels;
            for (size_t nIdx=0; nIdx < Pixel::NUM_NEIGHBORS; nIdx++)
            {
                Label * neighborLabel = neighborLabels[nIdx];
                // Only consider labels with valid values
                if (neighborLabel->value > 0)
                {
                    // Loop through uniqueLabels
                    bool unique = true;
                    for (size_t uIdx=0; uIdx < uniqueLabels.size(); uIdx++)
                    {
                        // Don't use the label if it's already been included
                        if (*neighborLabel == *(uniqueLabels[uIdx]))
                        {
                            unique = false;
                            break;
                        }
                    }
                    // If the label hasn't been added yet, do so
                    if (unique)
                    {
                        uniqueLabels.push_back(neighborLabel);
                    }
                }
            }
            
            // Determine what to do based on how many valid unique labels are nearby
            if (uniqueLabels.size() == 0)
            {
                //            std::cout << "tryToLabel: C.A" << std::endl;
                return false;
            }
            else if (uniqueLabels.size() == 1)
            {
                //            std::cout << "tryToLabel: C.B" << std::endl;
                label = uniqueLabels[0];
                label->addToPixelCount(1);
                return true;
            }
            else
            {
                //            std::cout << "tryToLabel: C.C" << std::endl;
                /**
                 For now, labels meeting the size criteria eat up the smaller ones.
                 This means that if there are 3 valid labels, but only one meets the size criteria,
                 then it will merge with one of the smaller labels before the 2 smaller labels
                 can merge with each other (and potentially create a valid sized region)
                 In the case of a tie between two large enough regions and one or more smaller, the merger will
                 occur arbitrarily (i.e., whichever bigger one is nearest in the list to the smaller merges with it)
                 **/
                
                // First loop: Big labels eat up smaller
                // Iterate from 1st idx so first iteration of loop doesn't cause problems
                bool prevBigEnough = uniqueLabels[0]->pixelCount >= Label::MIN_PIXEL_COUNT;
                for (size_t uIdx=1; uIdx < uniqueLabels.size(); uIdx++)
                {
                    // If the previous label does / doesn't meet the size criteria and this label doesn't / does, merge them
                    bool currBigEnough = uniqueLabels[uIdx]->pixelCount >= Label::MIN_PIXEL_COUNT;
                    if (prevBigEnough != currBigEnough)
                    {
                        uniqueLabels[uIdx]->mergeWithLabel(uniqueLabels[uIdx-1]);
                        std::cout << "Merging (A) : " << uniqueLabels[uIdx]->value << " with " << uniqueLabels[uIdx-1]->value << std::endl;
                    }
                    prevBigEnough = uniqueLabels[uIdx]->pixelCount >= Label::MIN_PIXEL_COUNT;
                }
                
                // Second loop: Smaller labels merge together
                prevBigEnough = uniqueLabels[0]->pixelCount >= Label::MIN_PIXEL_COUNT;
                for (size_t uIdx=1; uIdx < uniqueLabels.size(); uIdx++)
                {
                    // If the previous label doesn't meet the size criteria and this label doesn't, merge them
                    bool currBigEnough = uniqueLabels[uIdx]->pixelCount >= Label::MIN_PIXEL_COUNT;
                    if ( (!prevBigEnough) && (!currBigEnough) )
                    {
                        uniqueLabels[uIdx]->mergeWithLabel(uniqueLabels[uIdx-1]);
                        std::cout << "Merging (B) : " << uniqueLabels[uIdx]->value << " with " << uniqueLabels[uIdx-1]->value << std::endl;
                    }
                    prevBigEnough = uniqueLabels[uIdx]->pixelCount >= Label::MIN_PIXEL_COUNT;
                }
                
                return true;
            }
        }
    }
    
    bool isLabeled()
    {
        return label->value > 0;
    }
};

    /***************************************
     Watershed: Performs the 3D watershed
     **************************************/

class Watershed
{
public:
    // ATTRIBUTES
    std::vector<PixelType> * gridVals; // Input grid values (flattened)
    std::vector<PixelType> * seedVals;    // Input precomputed labels ([0,n]), no skips. 0 == no label, > 0 == label
    PixelType stepSize;                // Input step size used in creating levelSet
    PixelType lThresh;                 // Input lower threshold for levelSet
    PixelType uThresh;                 // Input upper threshold for levelSet
    int iSize;                      // Input size of grid in i dimension
    int jSize;                      // Input size of grid in j dimension
    int kSize;                      // Input size of grid in k dimension
    
    std::vector<Pixel> gridPixels;  // Grid pixels
    std::vector<Label> gridLabels;  // Container for Label objects used in grid
    LevelSet ls;                    // LevelSet object
    
    // METHODS
    /*
     Watershed(std::vector<double> *, std::vector<int> *, double, double, double, int, int, int);    // Constructor
     
     void setupInitialLabels();      // Setup initial labels and input label seeds
     void setupGridPixels();         // Create grid pixels
     void growWatershed();           // Grow entire watershed
     void growWatershedLevel(int);   // Grow watershed on single level of levelSet
     */
    
    Watershed(std::vector<PixelType> * wGridVals, std::vector<PixelType> * wSeedVals, double wStepSize, double wLThresh, double wUThresh, int wISize, int wJSize, int wKSize) :
    ls(wGridVals, wStepSize, wLThresh, wUThresh)
    {
        /**
         seedVals should range from [0,n], where 0 indicates an unlabeled pixel and n > 0 indicates a labeled pixel.  No gaps are allowed (e.g., [0,1,2,4] == broken)
         
         There are currently no true peaks in this implementation, as the seedVals make determining peaks from predetermined
         labels rather difficult, especially if some other method / data was used
         **/
        
        std::cout << " --- BEGIN WATERSHED --- " << std::endl;
        
        // Set attributes
        gridVals = wGridVals;
        seedVals = wSeedVals;
        stepSize = wStepSize;
        lThresh = wLThresh;
        uThresh = wUThresh;
        iSize = wISize;
        jSize = wJSize;
        kSize = wKSize;
        
        setupInitialLabels();
        setupGridPixels();
        growWatershed();
        
        std::cout << " --- END WATERSHED --- " << std::endl;
    }
    
    void setupInitialLabels()
    {
        std::cout << "Setting up initial labels" << std::endl;
        // Set up initial labels and input label seeds
        gridLabels.reserve(100); // Just a guess on size
        std::set<PixelType> uniqueLabelVals(seedVals->begin(), seedVals->end());
        for (std::set<PixelType>::iterator it=uniqueLabelVals.begin(); it != uniqueLabelVals.end(); it++)
        {
            PixelType seedVal = * it;
            Label newLabel(seedVal);
            gridLabels.push_back(newLabel);
        }
    }
    
    void setupGridPixels()
    {
        std::cout << "Creating grid pixels" << std::endl;
        // Create grid pixels
        gridPixels.reserve(gridVals->size());
        for (size_t idx=0; idx < gridVals->size(); idx++)
        {
            PixelType seedIdx = (*seedVals)[idx];
            Label * gridLabel = & (gridLabels[seedIdx]);
            Pixel newPixel(idx, iSize, jSize, kSize, gridLabel);
            gridPixels.push_back(newPixel);
        }
        
        std::cout << "Pointing pixels to neighbor labels" << std::endl;
        // Point grid pixels to neighbor labels
        Pixel::setupNeighborMod(iSize, jSize, kSize);
        for (size_t idx=0; idx < gridPixels.size(); idx++)
        {
            if (!gridPixels[idx].onEdge)
            {
                for (int nMIdx=0; nMIdx < Pixel::neighborMod.size(); nMIdx++)
                {
                    int nIdx = idx + Pixel::neighborMod[nMIdx];
                    gridPixels[idx].addNeighborLabel(gridPixels[nIdx].label);
                }
            }
        }
    }
    
    void growWatershed()
    {
        std::cout << "Growing watershed..." << std::endl;
        // Loop through levels in levelSet
        for (int lIdx=ls.nBins-1; lIdx >= 0 < ls.nBins; lIdx--)
        {
            growWatershedLevel(lIdx);
        }
    }
    
    void growWatershedLevel(int lIdx)
    {
        std::cout << "Processing level " << lIdx << std::endl;
        
        std::vector<int> levelIdcs = ls.binIdcs[lIdx];
        bool levelFinished = false;
        while (!levelFinished)
        {
            // Label what can be labeled without adding new peaks
            bool changeMade = true;
            while (changeMade)
            {
                std::cout << "iter begin" << std::endl;
                changeMade = false;
                for (int lPIdx = 0; lPIdx < levelIdcs.size(); lPIdx++)
                {
                    int pIdx = levelIdcs[lPIdx];
                    changeMade = ( (gridPixels[pIdx].tryToLabel()) || (changeMade) );
                }
                std::cout << "iter end" << std::endl;
            }
            std::cout << "Checking for new peaks" << std::endl;
            // Add new peak
            PixelType largestUnlabeledValue = static_cast<PixelType>(-1.0*DBL_MAX);
            int largestUnlabeledIdx = -1;
            for (int lPIdx=0; lPIdx < levelIdcs.size(); lPIdx++)
            {
                int pIdx = levelIdcs[lPIdx];
                if ( (gridPixels[pIdx].isLabeled()) && ((*gridVals)[pIdx] > largestUnlabeledValue) && (!gridPixels[pIdx].onEdge) )
                {
                    largestUnlabeledValue = (*gridVals)[pIdx];
                    largestUnlabeledIdx = pIdx;
                }
            }
            // If we found a new peak, create it and the new label
            if (largestUnlabeledIdx != -1)
            {
                Label newLabel(Label::numLabels);   // Create new label
                gridLabels.push_back(newLabel);     // Store it
                gridPixels[largestUnlabeledIdx].label = &(gridLabels[gridLabels.size()-1]); // Point peak pixel to new label location
            }
            // Otherwise no changes have been made and there are no new peaks available => then the level is complete
            else
            {
                levelFinished = true;
            }
        }
    }
    
};


// Declaration of static variables
int Label::numLabels = 0;
std::vector<int> Pixel::neighborMod;


namespace itk
{
    template< class T>
    class LakWatershedFilter : public ImageToImageFilter< T, T >
    {
        
    public:
        /** Standard class typedefs. */
        typedef LakWatershedFilter				 Self;
        typedef ImageToImageFilter< T, T >	 Superclass;
        typedef SmartPointer< Self >		 Pointer;
        typedef SmartPointer< const Self >	 ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(LakWatershedFilter, ImageToImageFilter);
        
        /** Parameter personalization */
        void SetStep(float lvl)
        {
            stepSize = lvl;
        }
        void SetLowerThesh(float lt)
        {
            lThresh =  lt;
        }
        void SetUpperThesh(float lt)
        {
            uThresh =  lt;
        }
        
        
    protected:
        LakWatershedFilter()
        {
            
        }
        
        // Does the real work.
        virtual void GenerateData()
        {
            typename T::ConstPointer input = this->GetInput();
            typename T::Pointer output = this->GetOutput();
            
            // Set size for the buffer
            iSize = input->GetLargestPossibleRegion().GetSize()[0];
            jSize = input->GetLargestPossibleRegion().GetSize()[1];
            kSize = input->GetLargestPossibleRegion().GetSize()[2];
            
            // flatten from itkImage to std::vector
            // faster alternate: gridVals.assign( input->GetBufferPointer() , input->GetBufferPointer() + len);
            gridVals.clear();
            seedVals.clear();
            
            ConstIteratorType  iteratorInput( input, input->GetLargestPossibleRegion()  );
            iteratorInput.GoToBegin();
            while ( !iteratorInput.IsAtEnd() )
            {
                seedVals.push_back( 0.0 );
                gridVals.push_back( iteratorInput.Get() );
                ++iteratorInput;
            }

            
            Watershed ws (&gridVals, &seedVals, stepSize, lThresh, uThresh, iSize, jSize, kSize);
            
            
            // allocate output
            output->SetRegions( input->GetLargestPossibleRegion() );
            output->SetSpacing( input->GetSpacing() );
            output->SetOrigin( input->GetOrigin() );
            output->Allocate();
            
            
            int i=0;
            IteratorType  iteratorOutput( output, output->GetLargestPossibleRegion() );
            iteratorOutput.GoToBegin();
            while ( !iteratorOutput.IsAtEnd() )
            {
                iteratorOutput.Set( gridVals[i++] );
                ++iteratorOutput;
            }
        }
        
        LakWatershedFilter(const Self &);  //purposely not implemented
        void operator=(const Self &);    //purposely not implemented
        
        
        std::vector< PixelType > gridVals;
        std::vector< PixelType > seedVals;
        int iSize, jSize, kSize;
        double stepSize, lThresh, uThresh;
        
        typedef itk::ImageRegionIterator< T > IteratorType;
        typedef itk::ImageRegionConstIterator< T > ConstIteratorType;
        
    };
    
} //namespace ITK



#endif /* defined(____Watershed__) */
