/**
  @class Indices

  A very simple class for putting an event with a given value in the correct histogram bin

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */


#ifndef INDICES_H
#define INDICES_H

#include <vector>

using namespace std;

class Indices
{
	public:
		Indices();
		Indices( const vector< unsigned int > & InputBinNumber, const vector< double > & InputMinima, const vector< double > & InputMaxima );
		~Indices();

		//Return the bin index corresponding to a particular value (or set of values)
		unsigned int GetIndex( const vector< double > & InputValues );

		//Return the index in each dimension
		vector< unsigned int > GetNDimensionalIndex( const vector< double > & InputValues );
		vector< unsigned int > GetNDimensionalIndex( const unsigned int InputIndex );

		//Return the bin central value in each dimension
		vector< double > GetCentralValues( const vector< unsigned int > & InputIndices );
		vector< double > GetCentralValues( const unsigned int InputIndex );

		//Definitions for the indices
		unsigned int GetBinNumber();
		unsigned int GetBinNumber( const unsigned int DimensionIndex );
		vector< double > GetMinima();
		vector< double > GetMaxima();

	protected:
		vector< double > minima, maxima, binWidths;
		vector< unsigned int > numberOfBins;
		unsigned int numberOfDimensions;

};

#endif
