/**
  @class UniformIndices

  A class for putting an event with a given value in the correct histogram bin, assuming the histogram bins have uniform widths

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */


#ifndef UNIFORM_INDICES_H
#define UNIFORM_INDICES_H

#include "IIndexCalculator.h"
#include <vector>

using namespace std;

class UniformIndices : public IIndexCalculator
{
	public:
		UniformIndices();
		UniformIndices( vector< unsigned int > InputBinNumbers, vector< double > InputMinima, vector< double > InputMaxima );
		virtual ~UniformIndices();

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

		vector< vector< double >* > binLowEdges;
		vector< double > minima, maxima, binWidths;
		vector< unsigned int > numberOfBins;
		unsigned int numberOfDimensions;
		vector< vector< double > > binValueSums, binValueNormalisations;
};

#endif
