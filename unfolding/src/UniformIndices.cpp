/**
  @class UniformIndices

  A very simple class for putting an event with a given value in the correct histogram bin, assuming the histogram bins have uniform widths

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "UniformIndices.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>

//Default constructor
UniformIndices::UniformIndices()
{
}

//Constructor for fixed bin widths
UniformIndices::UniformIndices( vector< unsigned int > InputBinNumbers, vector< double > InputMinima, vector< double > InputMaxima )
{
	numberOfBins = InputBinNumbers;

	//Check for stupidity
	if ( InputMinima.size() != InputMaxima.size() || InputMinima.size() != InputBinNumbers.size() )
	{
		cerr << "Different numbers of minima(" << InputMinima.size() << "), maxima(" << InputMaxima.size() << ") or bin numbers(" << InputBinNumbers.size() << ") are provided" << endl;
		exit(1);
	}
	else
	{
		numberOfDimensions = InputMinima.size();
	}

	//Loop over dimensions
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
	{
		//Check max and min are the right way round
		if ( InputMaxima[ dimensionIndex ] > InputMinima[ dimensionIndex ] )
		{
			minima.push_back( InputMinima[ dimensionIndex ] );
			maxima.push_back( InputMaxima[ dimensionIndex ] );
		}
		else
		{
			minima.push_back( InputMaxima[ dimensionIndex ] );
			maxima.push_back( InputMinima[ dimensionIndex ] );
		}

		//Calculate the bin width
		binWidths.push_back( ( maxima[ dimensionIndex ] - minima[ dimensionIndex ] ) / (double)numberOfBins[ dimensionIndex ] );
	}

	//Make the vectors to hold the bin central values
	for ( unsigned int dimensionIndex = 0; dimensionIndex < InputBinNumbers.size(); dimensionIndex++ )
	{
		vector< double > oneDimensionDistribution( numberOfBins[ dimensionIndex ] + 2, 0.0 );

		binValueSums.push_back( oneDimensionDistribution );
		binValueNormalisations.push_back( oneDimensionDistribution );

		//Make the bin low edges
		for ( unsigned int binIndex = 0; binIndex <= numberOfBins[ dimensionIndex ]; binIndex++ )
		{
			oneDimensionDistribution[ binIndex ] = minima[ dimensionIndex ] + ( binIndex * binWidths[ dimensionIndex ] );
		}
		//Remove the last entry - vector too long
		oneDimensionDistribution.pop_back();

		//Store the low edges - pointers so they can last outside this class instance
		binLowEdges.push_back( new vector< double >( oneDimensionDistribution ) );
	}
}

//Return a copy of the object, minus any stored data
IIndexCalculator * UniformIndices::Clone()
{
	return new UniformIndices( numberOfBins, minima, maxima );
}

//Destructor
UniformIndices::~UniformIndices()
{
	//Let the bin edge vectors exist past the end of the class
	//binLowEdges.clear();

	minima.clear();
	maxima.clear();
	binWidths.clear();
	numberOfBins.clear();
}

//Return the bin number / index corresponding to a given value
unsigned int UniformIndices::GetIndex( const vector< double > & InputValues )
{
	//Stupidity check
	if ( numberOfDimensions != InputValues.size() )
	{
		cerr << "Using a " << InputValues.size() << "D index lookup for a " << numberOfDimensions << "D unfolding" << endl;
		exit(1);
	}
	else
	{
		//Loop over all dimensions
		unsigned int totalIndex = 0;
		int binMultiplier = 1;
		for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
		{
			//Add to the total index
			totalIndex += GetOneDimensionIndex( InputValues[ dimensionIndex ], dimensionIndex ) * (double)binMultiplier;
			binMultiplier *= ( numberOfBins[ dimensionIndex ] + 2 );
		}

		return totalIndex;
	}
}

//Return the bin number in each dimension
vector< unsigned int > UniformIndices::GetNDimensionalIndex( const vector< double > & InputValues )
{
	//Stupidity check
	if ( numberOfDimensions != InputValues.size() )
	{
		cerr << "Using a " << InputValues.size() << "D index lookup for a " << numberOfDimensions << "D unfolding" << endl;
		exit(1);
	}
	else
	{
		//Loop over all dimensions
		vector< unsigned int > overallIndex( numberOfDimensions, 0 );
		for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
		{
			//Store the index
			overallIndex[ dimensionIndex ] = GetOneDimensionIndex( InputValues[ dimensionIndex ], dimensionIndex );
		}

		return overallIndex;
	}
}

//Return the bin number in each dimension
vector< unsigned int > UniformIndices::GetNDimensionalIndex( unsigned int InputIndex )
{
	//Loop over all dimensions
	vector< unsigned int > overallIndex( numberOfDimensions, 0 );
	unsigned int remainingIndexBits = InputIndex;
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
	{
		//Get the index in  this dimension
		unsigned int thisIndex = remainingIndexBits % ( numberOfBins[ dimensionIndex ] + 2 );
		overallIndex[ dimensionIndex ] = thisIndex;

		//Work out the remainder for the next dimension
		remainingIndexBits -= thisIndex;
		remainingIndexBits /= ( numberOfBins[ dimensionIndex ] + 2 );
	}

	return overallIndex;
}

//Return the central value of a given bin
vector< double > UniformIndices::GetCentralValues( const vector< unsigned int > & InputIndices )
{
	if ( numberOfDimensions != InputIndices.size() )
	{
		cerr << "Using a " << InputIndices.size() << "D index lookup for a " << numberOfDimensions << "D unfolding" << endl;
		exit(1);
	}
	else
	{
		//Loop over all dimensions
		vector< double > centralValues( numberOfDimensions, 0.0 );
		for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
		{
			unsigned int inputIndex = InputIndices[ dimensionIndex ];

			if ( inputIndex <= numberOfBins[ dimensionIndex ] + 1 )
			{
				//Find out the total of the values of all the events in this bin
				double binTotal = binValueSums[ dimensionIndex ][ inputIndex ];

				//Check if there's any data in the bin
				if ( binTotal == 0.0 )
				{
					//If the bin is empty, fall back to the old method
					if ( inputIndex == 0 )
					{
						//Underflow bin
						centralValues[ dimensionIndex ] =  minima[ dimensionIndex ] - ( 0.5 * binWidths[ dimensionIndex ] );
					}
					else if ( inputIndex == numberOfBins[ dimensionIndex ] + 1 )
					{
						//Overflow bin
						centralValues[ dimensionIndex ] = maxima[ dimensionIndex ] + ( 0.5 * binWidths[ dimensionIndex ] );
					}
					else if ( inputIndex <= numberOfBins[ dimensionIndex ] )
					{
						//Regular bin
						centralValues[ dimensionIndex ] = ( ( (double)inputIndex - 0.5 ) * binWidths[ dimensionIndex ] ) + minima[ dimensionIndex ];
					}
				}
				else
				{
					//Return the average value of the events in this bin
					centralValues[ dimensionIndex ] = binTotal / binValueNormalisations[ dimensionIndex ][ inputIndex ];
				}
			}
			else
			{
				cerr << "ERROR: Nonsensical index provided (" << inputIndex << ") when should be in range 0 to " << numberOfBins[ dimensionIndex ] + 1 << endl;
				exit(1);
			}
		}
		return centralValues;
	}
}
vector< double > UniformIndices::GetCentralValues( unsigned int InputIndex )
{
	//Work out the N-dim index first
	return GetCentralValues( GetNDimensionalIndex( InputIndex ) );
}

//Return the total bin number
unsigned int UniformIndices::GetBinNumber()
{
	unsigned int binNumber = 1;
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
	{
		binNumber *= ( numberOfBins[ dimensionIndex ] + 2 );
	}
	return binNumber;
}

//Return the number of bins in a given dimension
unsigned int UniformIndices::GetBinNumber( unsigned int DimensionIndex )
{
	if ( DimensionIndex >= numberOfDimensions )
	{
		cerr << "Specified invalid dimension index (" << DimensionIndex << ") - not in range 0 to " << numberOfDimensions - 1 << endl;
		exit(1);
	}
	else
	{
		return numberOfBins[ DimensionIndex ] + 2;
	}
}

//Return the index in a given dimension
unsigned int UniformIndices::GetOneDimensionIndex( double Value, unsigned int Dimension )
{
	if ( Value < minima[ Dimension ] )
	{
		//Underflow bin
		return 0;
	}
	else if ( Value > maxima[ Dimension ] )
	{
		//Overflow bin
		return numberOfBins[ Dimension ] + 1;
	}
	else
	{
		//Regular bin
		return ceil( ( Value - minima[ Dimension ] ) / binWidths[ Dimension ] );
	}
}

//Input the data to calculate the central values
void UniformIndices::StoreDataValue( const vector< double > & Data, double Weight )
{
	//Find out the bin that the event falls into in each dimension
	vector< unsigned int > binIndices = GetNDimensionalIndex( Data );

	//Loop over dimensions
	for ( unsigned int dimensionIndex = 0; dimensionIndex < binIndices.size(); dimensionIndex++ )
	{
		//Find out the index of the bin the event falls into in this dimension
		unsigned int binIndex = binIndices[dimensionIndex];

		//Add the event to the bin in this dimension
		binValueSums[ dimensionIndex ][ binIndex ] += Data[ dimensionIndex ] * Weight;
		binValueNormalisations[ dimensionIndex ][ binIndex ] += Weight;
	}
}

double * UniformIndices::GetBinLowEdgesForRoot( unsigned int DimensionIndex )
{
	return &( *( binLowEdges[ DimensionIndex ] ) )[0];
}
