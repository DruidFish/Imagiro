/**
  @class SmearingMatrix

  The crucial part of the unfolding process - the matrix describing detector effects on data

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "SmearingMatrix.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

//Default constructor
SmearingMatrix::SmearingMatrix()
{
}

//Constructor from MC data
SmearingMatrix::SmearingMatrix( IIndexCalculator * InputIndices )
{
	isFinalised = false;
	indexCalculator = InputIndices;

	//Initialise the normalisation
	normalisation = vector< double >( indexCalculator->GetBinNumber() + 1, 0.0 );
	totalPaired = 0.0;
	totalMissed = 0.0;
	totalFake = 0.0;
}

//Populate the matrix with values from events
void SmearingMatrix::StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight, double RecoWeight )
{
	//Look up the indices of the truth and reco
	unsigned int truthIndex = indexCalculator->GetIndex( Truth );
	unsigned int recoIndex = indexCalculator->GetIndex( Reco );

	//Increment values
	AddToEntry( truthIndex, recoIndex, RecoWeight );
	normalisation[ truthIndex ] += TruthWeight;
	totalPaired += TruthWeight;
}
void SmearingMatrix::StoreUnreconstructedTruth( vector< double > Truth, double Weight )
{
	//Look up the index of the truth value
	unsigned int truthIndex = indexCalculator->GetIndex( Truth );
	unsigned int recoIndex = indexCalculator->GetBinNumber();

	//Increment values
	AddToEntry( truthIndex, recoIndex, Weight );
	normalisation[ truthIndex ] += Weight;
	totalMissed += Weight;
}
void SmearingMatrix::StoreReconstructedFake( vector< double > Reco, double Weight )
{
	//Look up the index of the reco value
	unsigned int truthIndex = indexCalculator->GetBinNumber();
	unsigned int recoIndex = indexCalculator->GetIndex( Reco );

	//Increment values
	AddToEntry( truthIndex, recoIndex, Weight );
	normalisation[ truthIndex ] += Weight;
	totalFake += Weight;
}

//Do a bunch of extra calculations that aren't necessary if you just want the raw smearing matrix
void SmearingMatrix::Finalise()
{
	if ( !isFinalised )
	{
		unsigned int binNumber = indexCalculator->GetBinNumber() + 1;

		//Prepare storage for the efficiencies and effect probabilities
		efficiencies = vector< double >( binNumber, 0.0 );

		//Finalise each filled element
		map< pair< unsigned int, unsigned int >, double >::iterator matrixIterator;
		for ( matrixIterator = matrix.begin(); matrixIterator != matrix.end(); matrixIterator++ )
		{
			//Read the cause index for this entry
			unsigned int causeIndex = matrixIterator->first.first;

			//Normalise the entry
			double value = matrixIterator->second / normalisation[ causeIndex ];
			matrix[ matrixIterator->first ] = value;

			//Add to the efficiency of this cause
			efficiencies[ causeIndex ] += value;
		}

		VectorsFromMap( binNumber );
		isFinalised = true;
	}
}

//Destructor
SmearingMatrix::~SmearingMatrix()
{
	if ( isFinalised )
	{
		efficiencies.clear();
	}
	normalisation.clear();
}

//Return efficiency of a given cause
double SmearingMatrix::GetEfficiency( unsigned int CauseIndex )
{
	//Check if efficiencies are already calculated
	if ( !isFinalised )
	{
		Finalise();
	}

	return efficiencies[ CauseIndex ];
}

//Return total number of truth events in a given bin
double SmearingMatrix::GetTruthTotal( unsigned int CauseIndex )
{
	return normalisation[ CauseIndex ];
}

double SmearingMatrix::GetTotalPaired()
{
	return totalPaired;
}

double SmearingMatrix::GetTotalMissed()
{
	return totalMissed;
}

double SmearingMatrix::GetTotalFake()
{
	return totalFake;
}
