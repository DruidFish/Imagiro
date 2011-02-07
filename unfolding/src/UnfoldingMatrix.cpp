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
	//Calculate the probabilities of the effects
	int binNumber = InputIndices->GetBinNumber();
	vector<double> effectProbabilities;
	for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
	{
		double effectProbability = 0.0;

		for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
		{
			effectProbability += ( InputSmearing->GetElement( causeIndex, effectIndex ) * InputDistribution->GetBinProbability(causeIndex) );
		}

		effectProbabilities.push_back(effectProbability);
	}

	//Calculate the matrix elements
	matrix = vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) );
	for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		vector<double> causeElements;

		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			double numerator = InputSmearing->GetElement( causeIndex, effectIndex ) * InputDistribution->GetBinProbability(causeIndex);
			numerator /= ( effectProbabilities[effectIndex] * InputSmearing->GetEfficiency(causeIndex) );

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

	//Calculate the unfolding covariances
	/*covariance = vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) );
	for ( int indexU = 0; indexU < binNumber; indexU++ )
	{
		double efficiencyU = InputSmearing->GetEfficiency(indexU);

		for ( int indexI = 0; indexI < binNumber; indexI++ )
		{
			double justUIPart = matrix[indexU][indexI] * efficiencyU / InputSmearing->GetElement( indexU, indexI );

			//Check for divide by zero errors
			if ( isinf(justUIPart) )
			{
				cerr << "UNFOLDING: Divide by zero error (smearing)" << endl;
				exit(1);
			}
			else if ( isnan(justUIPart) )
			{
				justUIPart = 0.0;
			}

			covariance[indexI][indexU] = justUIPart;
		}
	}*/
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

//Return the unfolding covariance
double UnfoldingMatrix::GetCovarianceTerm( int IndexI, int IndexU )
{
	return covariance[IndexI][IndexU];
}
