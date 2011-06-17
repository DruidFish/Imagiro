/**
  @class Indices

  A class for putting an event with a given value in the correct histogram bin, if the histogram has custum bin widths

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-06-2011
 */


#ifndef CUSTOM_INDICES_H
#define CUSTOM_INDICES_H

#include "IIndexCalculator.h"
#include <vector>
#include <map>

using namespace std;

class CustomIndices : public IIndexCalculator
{
	public:
		CustomIndices();
		CustomIndices( const vector< vector< double > > & LowEdges );
		virtual ~CustomIndices();

		//Return a copy of the object, minus any stored data
                virtual IIndexCalculator * Clone();

		//Return the bin index corresponding to a particular value (or set of values)
		virtual unsigned int GetIndex( const vector< double > & InputValues );

		//Return the index in each dimension
		virtual vector< unsigned int > GetNDimensionalIndex( const vector< double > & InputValues );
		virtual vector< unsigned int > GetNDimensionalIndex( unsigned int InputIndex );

		//Return the bin central value in each dimension
		virtual vector< double > GetCentralValues( const vector< unsigned int > & InputIndices );
		virtual vector< double > GetCentralValues( unsigned int InputIndex );

		//Definitions for the indices
		virtual unsigned int GetBinNumber();
		virtual unsigned int GetBinNumber( unsigned int DimensionIndex );
		virtual double * GetBinLowEdgesForRoot( unsigned int DimensionIndex );

		//Input data to calculate the central values
                virtual void StoreDataValue( const vector< double > & Data, double Weight = 1.0 );

	private:
		unsigned int GetOneDimensionIndex( double Value, unsigned int Dimension );

		vector< vector< double > > binLowEdges;
		vector< vector< double >* > binLowEdgePointers;
		vector< unsigned int > numberOfBins;
		unsigned int numberOfDimensions;
		map< double, unsigned int > binHighEdgeMap;
		map< double, unsigned int >::iterator searchIterator;
		vector< vector< double > > binValueSums, binValueNormalisations;

};

#endif
