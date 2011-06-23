/**
  @class CustomIndices

  A class for putting an event with a given value in the correct histogram bin, if the histogram has custum bin widths

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-06-2011
 */

#include "CustomIndices.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>

using namespace std;

//Default constructor
CustomIndices::CustomIndices()
{
}

//Constructor for custom bin widths
CustomIndices::CustomIndices( const vector< vector< double > > & LowEdges )
{
	binLowEdges = LowEdges;
	numberOfDimensions = LowEdges.size();

	//Check the bins are in order
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
	{
		//Make a map for looking up the bin edges
		map< double, unsigned int > thisEdgeMap;

		//Check for adjacent bin edges being equal or not correctly ordered
		for ( unsigned int binIndex = 0; binIndex < LowEdges[ dimensionIndex ].size() - 1; binIndex++ )
		{
			double lowEdge = LowEdges[ dimensionIndex ][ binIndex ];
			double highEdge = LowEdges[ dimensionIndex ][ binIndex + 1 ];
			if ( lowEdge > highEdge )
			{
				cerr << "Input histogram bin edges are not ordered correctly: " << lowEdge << " > " << highEdge << endl;
				cerr << "To avoid ambiguity, please resolve this manually" << endl;
				exit(1);
			}
			else if ( lowEdge == highEdge )
			{
				cerr << "Two histogram bins given the same low edge: " << lowEdge << endl;
				cerr << "Either you copy-and-pasted incautiously, or you had a really clever idea you thought you'd try" << endl;
				cerr << "Regardless, it won't work: please fix this" << endl;
				exit(1);
			}
			else
			{
				thisEdgeMap[ highEdge ] = binIndex;
			}
		}

		//Store the map
		binHighEdgeMaps.push_back( thisEdgeMap );

		//Store the number of bins in this dimension (one fewer bin than edge)
		numberOfBins.push_back( LowEdges[ dimensionIndex ].size() - 1 );
	}

	//Make the vectors to hold the bin central values
	for ( unsigned int dimensionIndex = 0; dimensionIndex < LowEdges.size(); dimensionIndex++ )
	{
		//Make an vector of empty bins (include under and overflow)
		vector< double > oneDimensionDistribution( numberOfBins[ dimensionIndex ] + 2, 0.0 );

		binValueSums.push_back( oneDimensionDistribution );
		binValueNormalisations.push_back( oneDimensionDistribution );

		//Make pointers to vectors of bin low edges
		binLowEdgePointers.push_back( new vector< double >( LowEdges[ dimensionIndex ] ) );
	}
}

//Return a copy of the object, minus any stored data
IIndexCalculator * CustomIndices::Clone()
{
	return new CustomIndices( binLowEdges );
}

//Destructor
CustomIndices::~CustomIndices()
{
	//binLowEdges.clear();
	binHighEdgeMaps.clear();
	numberOfBins.clear();
}

//Return the bin number / index corresponding to a given value
unsigned int CustomIndices::GetIndex( const vector< double > & InputValues )
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
vector< unsigned int > CustomIndices::GetNDimensionalIndex( const vector< double > & InputValues )
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
vector< unsigned int > CustomIndices::GetNDimensionalIndex( unsigned int InputIndex )
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
vector< double > CustomIndices::GetCentralValues( const vector< unsigned int > & InputIndices )
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
						//Underflow bin, just pick lowest bin edge - 1
						centralValues[ dimensionIndex ] = binLowEdges[ dimensionIndex ][ 0 ] - 1.0;
					}
					else if ( inputIndex == numberOfBins[ dimensionIndex ] + 1 )
					{
						//Overflow bin, just pick highest bin edge + 1
						centralValues[ dimensionIndex ] = binLowEdges[ dimensionIndex ][ inputIndex - 1 ] + 1.0;
					}
					else
					{
						//Regular bin
						centralValues[ dimensionIndex ] = ( binLowEdges[ dimensionIndex ][ inputIndex - 1 ] + binLowEdges[ dimensionIndex ][ inputIndex ] ) * 0.5;
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
vector< double > CustomIndices::GetCentralValues( unsigned int InputIndex )
{
	//Work out the N-dim index first
	return GetCentralValues( GetNDimensionalIndex( InputIndex ) );
}

//Return the total bin number
unsigned int CustomIndices::GetBinNumber()
{
	unsigned int binNumber = 1;
	for ( unsigned int dimensionIndex = 0; dimensionIndex < numberOfDimensions; dimensionIndex++ )
	{
		binNumber *= ( numberOfBins[ dimensionIndex ] + 2 );
	}
	return binNumber;
}

//Return the number of bins in a given dimension
unsigned int CustomIndices::GetBinNumber( unsigned int DimensionIndex )
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
unsigned int CustomIndices::GetOneDimensionIndex( double Value, unsigned int Dimension )
{
	if ( Value <= binLowEdges[ Dimension ][ 0 ] )
	{
		//Underflow bin
		return 0;
	}
	else
	{
		//Returns an iterator pointing to the first element in the container whose key compares greater than the value
		searchIterator = binHighEdgeMaps[ Dimension ].upper_bound( Value );

		if ( searchIterator == binHighEdgeMaps[ Dimension ].end() )
		{
			//Overflow
			return numberOfBins[ Dimension ] + 1;
		}
		else
		{
			//Regular bin: add one beacuse of the underflow bin
			return searchIterator->second + 1;
		}
	}
}

//Input the data to calculate the central values
void CustomIndices::StoreDataValue( const vector< double > & Data, double Weight )
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

double * CustomIndices::GetBinLowEdgesForRoot( unsigned int DimensionIndex )
{
	return &( *( binLowEdgePointers[ DimensionIndex ] ) )[0];
}
