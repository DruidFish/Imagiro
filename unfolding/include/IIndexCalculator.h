/**
  @interface IIndexCalculator

  Classes for putting an event with a given value in the correct histogram bin

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-06-2011
 */


#ifndef I_INDEX_CALCULATOR_H
#define I_INDEX_CALCULATOR_H

#include <vector>

using namespace std;

class IIndexCalculator
{
	public:
		//Destructor
		virtual ~IIndexCalculator()
		{
		}

		//Return a copy of the object, minus any stored data
		virtual IIndexCalculator * Clone() = 0;

		//Return the bin index corresponding to a particular value (or set of values)
		virtual unsigned int GetIndex( const vector< double > & InputValues ) = 0;

		//Return the index in each dimension
		virtual vector< unsigned int > GetNDimensionalIndex( const vector< double > & InputValues ) = 0;
		virtual vector< unsigned int > GetNDimensionalIndex( unsigned int InputIndex ) = 0;

		//Return the bin central value in each dimension
		virtual vector< double > GetCentralValues( const vector< unsigned int > & InputIndices ) = 0;
		virtual vector< double > GetCentralValues( unsigned int InputIndex ) = 0;

		//Definitions for the indices
		virtual unsigned int GetBinNumber() = 0;
		virtual unsigned int GetBinNumber( unsigned int DimensionIndex ) = 0;
		virtual double * GetBinLowEdgesForRoot( unsigned int DimensionIndex ) = 0;

		//Input data to calculate the central values
                virtual void StoreDataValue( const vector< double > & Data, double Weight = 1.0 ) = 0;
};

#endif
