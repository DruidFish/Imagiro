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
		Indices( vector<int> InputBinNumber, vector<double> InputMinima, vector<double> InputMaxima );
		~Indices();

		//Return the bin index corresponding to a particular value (or set of values)
		int GetIndex( vector<double> InputValues );

		//Return the index in each dimension
		vector<int> GetNDimensionalIndex( vector<double> InputValues );
		vector<int> GetNDimensionalIndex( int InputIndex );

		//Return the bin central value in each dimension
		vector<double> GetCentralValues( vector<int> InputIndices );
		vector<double> GetCentralValues( int InputIndex );

		//Definitions for the indices
		int GetBinNumber();
		int GetBinNumber( int DimensionIndex );
		vector<double> GetMinima();
		vector<double> GetMaxima();

	protected:
		vector<double> minima, maxima, binWidths;
		vector<int> numberOfBins;

};

#endif
