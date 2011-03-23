/**
  @class DataIndices

  Extends the Indices class by using data to give more accurate bin central values  

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#include "DataIndices.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>

using namespace std;

//Default constructor
DataIndices::DataIndices()
{
}

//Constructor for N-dimensional indices
DataIndices::DataIndices( vector<int> InputBinNumbers, vector<double> InputMinima, vector<double> InputMaxima ) : Indices( InputBinNumbers, InputMinima, InputMaxima )
{
	//Make the vectors to hold the bin central values
	for ( unsigned int dimensionIndex = 0; dimensionIndex < InputBinNumbers.size(); dimensionIndex++ )
	{
		vector<double> oneDimensionDistribution( GetBinNumber( dimensionIndex ), 0.0 );

		binValueSums.push_back( oneDimensionDistribution );
		binValueNormalisations.push_back( oneDimensionDistribution );
	}
}

//Destructor
DataIndices::~DataIndices()
{
}

//Return the central value of a given bin
vector<double> DataIndices::GetCentralValuesFromData( vector<int> InputIndices )
{
	vector<double> centralValues( InputIndices.size(), 0.0 );

	//Loop over dimensions
	for ( unsigned int dimensionIndex = 0; dimensionIndex < InputIndices.size(); dimensionIndex++ )
	{
		//Find out the index of the bin we're returning in this dimension
		int binIndex = InputIndices[ dimensionIndex ];

		if ( binIndex >= 0 && binIndex <= numberOfBins[ dimensionIndex ] + 1)
		{
			//Find out the total of the values of all the events in this bin
			double binTotal = binValueSums[ dimensionIndex ][ binIndex ];

			//If the bin is empty, fall back to the old method
			if ( binTotal == 0.0 )
			{
				if ( binIndex == 0 )
				{
					//Underflow bin
					centralValues[ dimensionIndex ] = minima[ dimensionIndex ] - ( 0.5 * binWidths[ dimensionIndex ] );
				}
				else if ( binIndex == numberOfBins[ dimensionIndex ] + 1 )
				{
					//Overflow bin
					centralValues[ dimensionIndex ] = maxima[ dimensionIndex ] + ( 0.5 * binWidths[ dimensionIndex ] );
				}
				else
				{
					//Regular bin
					centralValues[ dimensionIndex ] = ( ( (double)binIndex - 0.5 ) * binWidths[ dimensionIndex ] ) + minima[ dimensionIndex ];
				}
			}
			else
			{
				//Return the average value of the events in this bin
				centralValues[ dimensionIndex ] = binTotal / binValueNormalisations[ dimensionIndex ][ binIndex ];
			}
		}
		else
		{
			cerr << "ERROR: Nonsensical index provided (" << binIndex << ") when should be in range 0 to " << numberOfBins[ dimensionIndex ] + 1 << endl;
			exit(1);
		}
	}

	return centralValues;
}
vector<double> DataIndices::GetCentralValuesFromData( int InputIndex )
{
	return GetCentralValuesFromData( GetNDimensionalIndex( InputIndex ) );
}

//Input the data to calculate the central values
void DataIndices::StoreDataValue( vector<double> Data, double Weight )
{
	//Find out the bin that the event falls into in each dimension
	vector<int> binIndices = GetNDimensionalIndex( Data );

	//Loop over dimensions
	for ( unsigned int dimensionIndex = 0; dimensionIndex < binIndices.size(); dimensionIndex++ )
	{
		//Find out the index of the bin the event falls into in this dimension
		int binIndex = binIndices[dimensionIndex];

		//Add the event to the bin in this dimension
		binValueSums[ dimensionIndex ][ binIndex ] += Data[ dimensionIndex ] * Weight;
		binValueNormalisations[ dimensionIndex ][ binIndex ] += Weight;
	}
}
