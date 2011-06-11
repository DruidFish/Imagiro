/**
  @class DataIndices

  Extends the Indices class by using data to give more accurate bin central values

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#ifndef DATA_INDICES_H
#define DATA_INDICES_H

#include "Indices.h"

using namespace std;

class DataIndices : public Indices
{
	public:
		DataIndices();
		DataIndices( vector< unsigned int > InputBinNumber, vector< double > InputMinima, vector< double > InputMaxima );
		~DataIndices();

		//Return the bin central value in each dimension
		vector< double > GetCentralValuesFromData( vector< unsigned int > InputIndices );
		vector< double > GetCentralValuesFromData( unsigned int InputIndex );

		//Input the data to calculate the central values
		void StoreDataValue( vector< double > Data, double Weight = 1.0 );

	private:
		vector< vector< double > > binValueSums, binValueNormalisations;

};

#endif
