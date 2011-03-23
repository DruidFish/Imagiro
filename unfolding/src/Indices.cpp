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
Indices::Indices( vector<int> InputBinNumber, vector<double> InputMinima, vector<double> InputMaxima ) : numberOfBins(InputBinNumber)
{
	//Check for stupidity
	if ( InputMinima.size() != InputMaxima.size() || InputMaxima.size() != InputBinNumber.size() )
	{
		cerr << "Different numbers of minima(" << InputMinima.size() << "), maxima(" << InputMaxima.size() << ") or bin numbers(" << InputBinNumber.size() << ") are provided" << endl;
		exit(1);
	}

	//Loop over dimensions
	for ( unsigned int dimensionIndex = 0; dimensionIndex < InputMinima.size(); dimensionIndex++ )
	{
		//Check for sensible bin number
		if ( numberOfBins[dimensionIndex] <= 0 )
		{
			cerr << "Impossible number of bins (" << numberOfBins[dimensionIndex] << ")" << endl;
			exit(1);
		}

		//Check max and min are the right way round
		if ( InputMaxima[dimensionIndex] > InputMinima[dimensionIndex] )
		{
			minima.push_back( InputMinima[dimensionIndex] );
			maxima.push_back( InputMaxima[dimensionIndex] );
		}
		else
		{
			minima.push_back( InputMaxima[dimensionIndex] );
			maxima.push_back( InputMinima[dimensionIndex] );
		}

		//Calculate the bin width
		binWidths.push_back( ( maxima[dimensionIndex] - minima[dimensionIndex] ) / (double)numberOfBins[dimensionIndex] );
	}
}

//Destructor
Indices::~Indices()
{
}

//Return the bin number / index corresponding to a given value
int Indices::GetIndex( vector<double> InputValues )
{
	if ( minima.size() != InputValues.size() )
	{
		cerr << "Using a " << InputValues.size() << "D index lookup for a " << minima.size() << "D unfolding" << endl;
		exit(1);
	}
	else
	{
		//Loop over all dimensions
		int totalIndex = 0;
		int binMultiplier = 1;
		for ( unsigned int dimensionIndex = 0; dimensionIndex < minima.size(); dimensionIndex++ )
		{
			int thisIndex;
			if ( InputValues[dimensionIndex] < minima[dimensionIndex] )
			{
				//Underflow bin
				thisIndex = 0;
			}
			else if ( InputValues[dimensionIndex] > maxima[dimensionIndex] )
			{
				//Overflow bin
				thisIndex = numberOfBins[dimensionIndex] + 1;
			}
			else
			{
				thisIndex = ceil( ( InputValues[dimensionIndex] - minima[dimensionIndex] ) / binWidths[dimensionIndex] );
			}

			totalIndex += thisIndex * (double)binMultiplier;
			binMultiplier *= ( numberOfBins[dimensionIndex] + 2 );
		}

		return totalIndex;
	}
}

//Return the bin number in each dimension
vector<int> Indices::GetNDimensionalIndex( vector<double> InputValues )
{
	if ( minima.size() != InputValues.size() )
	{
		cerr << "Using a " << InputValues.size() << "D index lookup for a " << minima.size() << "D unfolding" << endl;
		exit(1);
	}
	else
	{
		//Loop over all dimensions
		vector<int> overallIndex( minima.size(), 0 );
		for ( unsigned int dimensionIndex = 0; dimensionIndex < minima.size(); dimensionIndex++ )
		{
			int thisIndex;
			if ( InputValues[dimensionIndex] < minima[dimensionIndex] )
			{
				//Underflow bin
				thisIndex = 0;
			}
			else if ( InputValues[dimensionIndex] > maxima[dimensionIndex] )
			{
				//Overflow bin
				thisIndex = numberOfBins[dimensionIndex] + 1;
			}
			else
			{
				thisIndex = ceil( ( InputValues[dimensionIndex] - minima[dimensionIndex] ) / binWidths[dimensionIndex] );
			}

			overallIndex[ dimensionIndex ] = thisIndex;
		}

		return overallIndex;
	}
}

//Return the bin number in each dimension
vector<int> Indices::GetNDimensionalIndex( int InputIndex )
{
	//Loop over all dimensions
	vector<int> overallIndex( minima.size(), 0 );
	int remainingIndexBits = InputIndex;
	for ( unsigned int dimensionIndex = 0; dimensionIndex < minima.size(); dimensionIndex++ )
	{
		int thisIndex = remainingIndexBits % ( numberOfBins[ dimensionIndex ] + 2 );
		overallIndex[ dimensionIndex ] = thisIndex;
		remainingIndexBits -= thisIndex;
		remainingIndexBits /= ( numberOfBins[ dimensionIndex ] + 2 );
	}

	return overallIndex;
}

//Return the central value of a given bin
vector<double> Indices::GetCentralValues( vector<int> InputIndices )
{
	if ( minima.size() != InputIndices.size() )
	{
		cerr << "Using a " << InputIndices.size() << "D index lookup for a " << minima.size() << "D unfolding" << endl;
		exit(1);
	}
	else
	{
		vector<double> centralValues( minima.size(), 0.0 );
		for ( unsigned int dimension = 0; dimension < minima.size(); dimension++ )
		{
			int inputIndex = InputIndices[ dimension ];

			if ( inputIndex == 0 )
			{
				//Underflow bin
				centralValues[ dimension ] =  minima[dimension] - ( 0.5 * binWidths[dimension] );
			}
			else if ( inputIndex == numberOfBins[dimension] + 1 )
			{
				//Overflow bin
				centralValues[ dimension ] = maxima[dimension] + ( 0.5 * binWidths[dimension] );
			}
			else if ( inputIndex > 0 && inputIndex <= numberOfBins[dimension] )
			{
				centralValues[ dimension ] = ( ( (double)inputIndex - 0.5 ) * binWidths[dimension] ) + minima[dimension];
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
vector<double> Indices::GetCentralValues( int InputIndex )
{
	return GetCentralValues( GetNDimensionalIndex( InputIndex ) );
}

//Return the defining values
vector<double> Indices::GetMinima()
{
	return minima;
}
vector<double> Indices::GetMaxima()
{
	return maxima;
}
int Indices::GetBinNumber()
{
	int binNumber = 1;
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfBins.size(); dimensionIndex++ )
	{
		binNumber *= ( numberOfBins[dimensionIndex] + 2 );
	}
	return binNumber;
}
int Indices::GetBinNumber( unsigned int DimensionIndex )
{
	if ( DimensionIndex < 0 || DimensionIndex >= numberOfBins.size() )
	{
		cerr << "Specified invalid dimension index (" << DimensionIndex << ") - not in range 0 to " << numberOfBins.size() - 1 << endl;
		exit(1);
	}
	else
	{
		return numberOfBins[DimensionIndex] + 2;
	}
}
