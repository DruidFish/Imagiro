/**
  @class CovarianceMatrix

  An attempt to implement an efficient version of D'Agostini's error propagation using sparse matrices

  @author Benjamin M Wynne bwynne@cern.ch
  @date 08-04-2011
 */

#include "CovarianceMatrix.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

const int PERCENT_PROGRESS_INCREMENT = 10;

//Default constructor - useless
CovarianceMatrix::CovarianceMatrix()
{
}

//Constructor that will actually calculate the covariance matrix
CovarianceMatrix::CovarianceMatrix( UnfoldingMatrix * InputUnfolding, SmearingMatrix * InputSmearing, Distribution * DataDistribution, double CorrectedSum, bool JustVariance )
{
	//Store some input
	inputSmearing = InputSmearing;
	inputUnfolding = InputUnfolding;
	correctedSum = CorrectedSum;

	//Status message
	time(&timeNow);
	if ( JustVariance )
	{
		cout << endl << "Started calculating bin variance values: " << ctime( &timeNow ) << endl;
	}
	else
	{
		cout << endl << "Started calculating full covariance matrix: " << ctime( &timeNow ) << endl;
	}

	//Construct a sparse matrix for the smearing covariance
	rsuMatrix = new SmearingCovariance( InputSmearing, InputUnfolding );

	//The covariance matrix will have zero entries unless there is a product of two non-zero unfolding matrix entries
	//Apologies for lack of better index labels, but probably best just to stick with the notation in D'Agostini's paper
	unsigned int firstEntryNumber = InputUnfolding->GetEntryNumberAndResetIterator();
	unsigned int percentIncrement = ceil( (double)firstEntryNumber / (double)PERCENT_PROGRESS_INCREMENT );
	unsigned int currentPercentTarget = PERCENT_PROGRESS_INCREMENT;
	unsigned int currentEntryTarget = percentIncrement;
	cout << "PROGRESS:";
	for ( unsigned int firstEntryIndex = 0; firstEntryIndex < firstEntryNumber; firstEntryIndex++ )
	{
		///Get a non-zero unfolding matrix entry
		unsigned int k, i;
		double firstEntryValue = InputUnfolding->GetNextEntry( k, i );

		//Check that the corresponding data bin is non-zero
		double dataI = DataDistribution->GetBinNumber( i );
		if ( dataI > 0.0 )
		{
			//Multiply by the data value
			firstEntryValue *= dataI;

			//Check to see if we only need the diagonal elements
			if ( JustVariance )
			{
				//Get all entries with the same cause index
				vector< unsigned int > * jIndices = InputUnfolding->GetIndicesWithFirstIndex( k );
				vector< double > * secondEntryValues = InputUnfolding->GetEntriesWithFirstIndex( k );

				//Loop over these entries
				for ( unsigned int secondEntryIndex = 0; secondEntryIndex < jIndices->size(); secondEntryIndex++ )
				{
					unsigned int j = ( *jIndices )[ secondEntryIndex ];

					//Check that the corresponding data bin is non-zero
					double dataJ = DataDistribution->GetBinNumber( j );
					if ( dataJ > 0.0 )
					{
						//Do the calculation
						CovarianceCalculation( i, j, k, k, firstEntryValue * ( *secondEntryValues )[ secondEntryIndex ], dataI, dataJ );

					}
				}
			}
			else
			{
				//Loop over all non-zero entries
				unsigned int secondEntryNumber = InputUnfolding->GetEntryNumberAndResetIterator( true );
				for ( unsigned int secondEntryIndex = 0; secondEntryIndex < secondEntryNumber; secondEntryIndex++ )
				{
					unsigned int l, j;
					double secondEntryValue = InputUnfolding->GetNextEntry( l, j, true );

					//Check that the corresponding data bin is non-zero
					double dataJ = DataDistribution->GetBinNumber( j );
					if ( dataJ > 0.0 )
					{
						//Do the calculation
						CovarianceCalculation( i, j, k, l, firstEntryValue * secondEntryValue, dataI, dataJ );
					}
				}
			}
		}

		//Progress indicator
		if ( firstEntryIndex == currentEntryTarget )
		{
			cout << " " << currentPercentTarget << "%";
			cout.flush();
			currentPercentTarget += PERCENT_PROGRESS_INCREMENT;
			currentEntryTarget += percentIncrement;
		}
	}
	VectorsFromMap( InputUnfolding->GetBinNumber() );

	//Status message
	time(&timeNow);
	cout << endl << endl << "Finished calculating covariance: " << ctime( &timeNow ) << endl;
}

void CovarianceMatrix::CovarianceCalculation( unsigned int I, unsigned int J, unsigned int K, unsigned int L, double unfoldingProductTimesDataI, double dataI, double dataJ )
{
	double covarianceContribution = 0.0;
	double mostOfIt = unfoldingProductTimesDataI * dataJ;

	//Smearing error contribution
	covarianceContribution += mostOfIt * rsuMatrix->ThisContribution( I, J, K, L );

	//Data error contribution
	if ( I == J )
	{
		covarianceContribution += unfoldingProductTimesDataI * ( 1.0 - ( dataI / correctedSum ) );
	}
	else
	{
		covarianceContribution -= mostOfIt / correctedSum;
	}

	//Store the result
	AddToEntry( K, L, covarianceContribution );
}

CovarianceMatrix::~CovarianceMatrix()
{
}
