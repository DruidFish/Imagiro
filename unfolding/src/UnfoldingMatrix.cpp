/**
  @class UnfoldingMatrix

  The Bayesian "inverse" of the smearing matrix, used to unfold the data
  Sparse version thereof

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "UnfoldingMatrix.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

//Default constructor
UnfoldingMatrix::UnfoldingMatrix()
{
}

//Constructor with correct arguments
UnfoldingMatrix::UnfoldingMatrix( SmearingMatrix * InputSmearing, Distribution * InputDistribution )
{
	//Get the dimension of the matrix
	int binNumber = InputSmearing->GetBinNumber();

	//Calculate the probabilities of the effects
	vector<double> effectProbabilities( binNumber, 0.0 );
	int entryNumber = InputSmearing->GetEntryNumberAndResetIterator();
	for ( int smearingIndex = 0; smearingIndex < entryNumber; smearingIndex++ )
	{
		//Retrieve the smearing matrix entry
		int causeIndex, effectIndex;
		double smearingValue = InputSmearing->GetNextEntry( causeIndex, effectIndex );

		//Add to the effect probability
		effectProbabilities[ effectIndex ] += smearingValue * InputDistribution->GetBinProbability( causeIndex );
	}

	//Calculate the matrix elements
	entryNumber = InputSmearing->GetEntryNumberAndResetIterator();
	for ( int smearingIndex = 0; smearingIndex < entryNumber; smearingIndex++ )
	{
		//Retrieve the smearing matrix entry    
		int causeIndex, effectIndex;            
		double smearingValue = InputSmearing->GetNextEntry( causeIndex, effectIndex );

		//Ignore zero entries
		double causeProbability = InputDistribution->GetBinProbability( causeIndex );
		if ( causeProbability != 0.0 )
		{
			//Work out the unfolding matrix entry
			double numerator = smearingValue * causeProbability;
			numerator /= ( effectProbabilities[ effectIndex ] * InputSmearing->GetEfficiency( causeIndex ) );

			//Store the entry
			AddToEntry( causeIndex, effectIndex, numerator );
		}
	}

	VectorsFromMap( binNumber );
}

//Destructor
UnfoldingMatrix::~UnfoldingMatrix()
{
}
