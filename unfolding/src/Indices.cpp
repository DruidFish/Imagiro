/**
  @class Indices

  A very simple class for putting an event with a given value in the correct histogram bin

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "Indices.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>

using namespace std;

//Default constructor
Indices::Indices()
{
}

//Constructor for N-dimensional indices
Indices::Indices( const vector< unsigned int > & InputBinNumber, const vector< double > & InputMinima, const vector< double > & InputMaxima )
{
	numberOfBins = InputBinNumber;

	//Check for stupidity
	if ( InputMinima.size() != InputMaxima.size() || InputMinima.size() != InputBinNumber.size() )
	{
		cerr << "Different numbers of minima(" << InputMinima.size() << "), maxima(" << InputMaxima.size() << ") or bin numbers(" << InputBinNumber.size() << ") are provided" << endl;
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
}

//Destructor
Indices::~Indices()
{
	minima.clear();
	maxima.clear();
	binWidths.clear();
	numberOfBins.clear();
}

//Return the bin number / index corresponding to a given value
unsigned int Indices::GetIndex( const vector< double > & InputValues )
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
			unsigned int thisIndex;
			if ( InputValues[ dimensionIndex ] < minima[ dimensionIndex ] )
			{
				//Underflow bin
				thisIndex = 0;
			}
			else if ( InputValues[ dimensionIndex ] > maxima[ dimensionIndex ] )
			{
				//Overflow bin
				thisIndex = numberOfBins[ dimensionIndex ] + 1;
			}
			else
			{
				//Regular bin
				thisIndex = ceil( ( InputValues[ dimensionIndex ] - minima[ dimensionIndex ] ) / binWidths[ dimensionIndex ] );
			}

			//Add to the total index
			totalIndex += thisIndex * (double)binMultiplier;
			binMultiplier *= ( numberOfBins[ dimensionIndex ] + 2 );
		}

		return totalIndex;
	}
}

//Return the bin number in each dimension
vector< unsigned int > Indices::GetNDimensionalIndex( const vector< double > & InputValues )
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
			unsigned int thisIndex;
			if ( InputValues[ dimensionIndex ] < minima[ dimensionIndex ] )
			{
				//Underflow bin
				thisIndex = 0;
			}
			else if ( InputValues[ dimensionIndex ] > maxima[ dimensionIndex ] )
			{
				//Overflow bin
				thisIndex = numberOfBins[ dimensionIndex ] + 1;
			}
			else
			{
				//Regular bin
				thisIndex = ceil( ( InputValues[ dimensionIndex ] - minima[ dimensionIndex ] ) / binWidths[ dimensionIndex ] );
			}

			//Store the index
			overallIndex[ dimensionIndex ] = thisIndex;
		}

		return overallIndex;
	}
}

//Return the bin number in each dimension
vector< unsigned int > Indices::GetNDimensionalIndex( const unsigned int InputIndex )
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
vector< double > Indices::GetCentralValues( const vector< unsigned int > & InputIndices )
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
		for ( unsigned int dimension = 0; dimension < numberOfDimensions; dimension++ )
		{
			unsigned int inputIndex = InputIndices[ dimension ];

			if ( inputIndex == 0 )
			{
				//Underflow bin
				centralValues[ dimension ] =  minima[ dimension ] - ( 0.5 * binWidths[ dimension ] );
			}
			else if ( inputIndex == numberOfBins[ dimension ] + 1 )
			{
				//Overflow bin
				centralValues[ dimension ] = maxima[ dimension ] + ( 0.5 * binWidths[ dimension ] );
			}
			else if ( inputIndex <= numberOfBins[ dimension ] )
			{
				//Regular bin
				centralValues[ dimension ] = ( ( (double)inputIndex - 0.5 ) * binWidths[ dimension ] ) + minima[ dimension ];
			}
			else
			{
				cerr << "ERROR: Nonsensical index provided (" << inputIndex << ") when should be in range 0 to " << numberOfBins[ dimension ] + 1 << endl;
				exit(1);
			}
		}
		return centralValues;
	}
}
vector< double > Indices::GetCentralValues( unsigned int InputIndex )
{
	//Work out the N-dim index first
	return GetCentralValues( GetNDimensionalIndex( InputIndex ) );
}

//Return the defining values
vector< double > Indices::GetMinima()
{
	return minima;
}
vector< double > Indices::GetMaxima()
{
	return maxima;
}
unsigned int Indices::GetBinNumber()
{
	unsigned int binNumber = 1;
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
	{
		binNumber *= ( numberOfBins[ dimensionIndex ] + 2 );
	}
	return binNumber;
}
unsigned int Indices::GetBinNumber( unsigned int DimensionIndex )
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
