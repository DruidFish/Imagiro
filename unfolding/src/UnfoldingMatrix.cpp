/**
  @class UnfoldingMatrix

  The Bayesian "inverse" of the smearing matrix, used to unfold the data

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
UnfoldingMatrix::UnfoldingMatrix( SmearingMatrix * InputSmearing, Distribution * InputDistribution, Indices * InputIndices )
{
	//Get the dimension of the matrix (include a bad bin)
	int binNumber = InputIndices->GetBinNumber() + 1;

	//Keep count of the number of events in the source distribution
	double totalCauseProbability = 0.0;

	//Calculate the probabilities of the effects
	vector<double> effectProbabilities( binNumber, 0.0 );
	for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		//Find out the probability of a particular cause
		double causeProbability = InputDistribution->GetBinProbability( causeIndex );

		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			//Add to the effect probability
			effectProbabilities[ effectIndex ] += ( InputSmearing->GetElement( causeIndex, effectIndex ) * causeProbability );
		}
	}

	//Calculate the matrix elements
	matrix = vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) );
	for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		//Find out the probability of a particular cause
		double causeProbability = InputDistribution->GetBinProbability( causeIndex );

		//Work out the corresponding matrix elements
		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			double numerator = InputSmearing->GetElement( causeIndex, effectIndex ) * causeProbability;
			numerator /= ( effectProbabilities[ effectIndex ] * InputSmearing->GetEfficiency( causeIndex ) );

			//Check for divide-by-zero issues
			if ( isinf(numerator) )
			{
				cerr << "UNFOLDING: Divide by zero error (cause/effect)" << endl;
				exit(1);
			}
			else if ( isnan(numerator) )
			{
				numerator = 0.0;
			}

			matrix[causeIndex][effectIndex] = numerator;
		}
	}
}

//Destructor
UnfoldingMatrix::~UnfoldingMatrix()
{
	matrix.clear();
}

//Return the matrix element
double UnfoldingMatrix::GetElement( int CauseIndex, int EffectIndex )
{
	return matrix[CauseIndex][EffectIndex];
}
