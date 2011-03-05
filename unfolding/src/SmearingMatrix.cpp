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
SmearingMatrix::SmearingMatrix( Indices * InputIndices ) : indexCalculator(InputIndices), isFinalised(false), totalPaired(0.0),
	totalMissed(0.0), totalFake(0.0)
{
	//Initialise the matrix and normalisation
	int binNumber = indexCalculator->GetBinNumber() + 1;
	matrix = vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) );
	normalisation = vector<double>( binNumber, 0.0 );
}

//Populate the matrix with unnormalised input smearing matrices to combine
void SmearingMatrix::StoreUnnormalisedMatrix( TH2F * InputMatrix )
{
	cerr << "This won't work" << endl;
	exit(1);

	//Check the matrix is the right size
	int binNumber = indexCalculator->GetBinNumber() + 1;
	int xBins = InputMatrix->GetNbinsX();
	int yBins = InputMatrix->GetNbinsY();
	if ( xBins != binNumber || yBins != binNumber)
	{
		cerr << "Input smearing matrix has incorrect dimensions: " << xBins << ", " << yBins << " vs " << binNumber << ", " << binNumber << endl;
		exit(1);
	}

	//Populate the matrix and normalisation
	for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			//Find the bin number
			int inputBinNumber = InputMatrix->GetBin( causeIndex + 1, effectIndex + 1, 0 );

			//Retrieve the value in the bin
			double inputBinValue = InputMatrix->GetBinContent(inputBinNumber);

			//Store the value
			matrix[causeIndex][effectIndex] += inputBinValue;
			normalisation[causeIndex] += inputBinValue;
		}
	}
}

//Populate the matrix with values from events
void SmearingMatrix::StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double TruthWeight, double RecoWeight )
{
	//Look up the indices of the truth and reco
	int truthIndex = indexCalculator->GetIndex(Truth);
	int recoIndex = indexCalculator->GetIndex(Reco);

	//Increment values
	matrix[truthIndex][recoIndex] += RecoWeight;
	normalisation[truthIndex] += TruthWeight;
	totalPaired += TruthWeight;
}
void SmearingMatrix::StoreUnreconstructedTruth( vector<double> Truth, double Weight )
{
	//Look up the index of the truth value
	int truthIndex = indexCalculator->GetIndex(Truth);
	int recoIndex = indexCalculator->GetBinNumber();

	//Increment values
	matrix[truthIndex][recoIndex] += Weight;
	normalisation[truthIndex] += Weight;
	totalMissed += Weight;
}
void SmearingMatrix::StoreReconstructedFake( vector<double> Reco, double Weight )
{
	//Look up the index of the reco value
	int truthIndex = indexCalculator->GetBinNumber();
	int recoIndex = indexCalculator->GetIndex(Reco);

	//Increment values
	matrix[truthIndex][recoIndex] += Weight;
	normalisation[truthIndex] += Weight;
	totalFake += Weight;
}

//Do a bunch of extra calculations that aren't necessary if you just want the raw smearing matrix
void SmearingMatrix::Finalise()
{
	if ( !isFinalised )
	{
		//Find matrix dimensions (include a bad bin)
		int binNumber = indexCalculator->GetBinNumber() + 1;

		//Prepare storage for the efficiencies
		efficiencies = vector<double>( binNumber, 0.0 );

		for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
		{
			//If there's no information, just make a uniform row
			if ( normalisation[ causeIndex ] == 0.0 )
			{
				for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
				{
					matrix[ causeIndex ][ effectIndex ] = matrix[ binNumber - 1 ][ effectIndex ] / (double)( binNumber - 1 );
				}

				efficiencies[ causeIndex ] = normalisation[ binNumber - 1 ] / (double)( binNumber - 1 );
			}
			else
			{
				for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
				{
					//Add the fakes
					if ( causeIndex != binNumber - 1 )
					{
						matrix[ causeIndex ][ effectIndex ] += ( matrix[ binNumber - 1 ][ effectIndex ] / (double)( binNumber - 1 ) );
					}

					//Calculate the efficiency
					efficiencies[ causeIndex ] += matrix[ causeIndex ][ effectIndex ];

					//Normalise
					matrix[ causeIndex ][ effectIndex ] /= normalisation[ causeIndex ];
				}

				efficiencies[ causeIndex ] /= normalisation[ causeIndex ];
			}
		}

		isFinalised = true;
	}
}

//Destructor
SmearingMatrix::~SmearingMatrix()
{
	matrix.clear();
}

//Return matrix element by index
double SmearingMatrix::GetElement( int CauseIndex, int EffectIndex )
{
	//Check if matrix already normalised
	if ( !isFinalised )
	{
		Finalise();
	}

	return matrix[ CauseIndex ][ EffectIndex ];
}

//Return efficiency of a given cause
double SmearingMatrix::GetEfficiency( int CauseIndex )
{
	//Check if efficiencies are already calculated
	if ( !isFinalised )
	{
		Finalise();
	}

	return efficiencies[CauseIndex];
}

//Return ratio of fake events to paired + missed
double SmearingMatrix::GetFakeRatio()
{
	return totalFake / ( totalPaired + totalMissed );
}

//Return ratio of missed events to paired + fake
double SmearingMatrix::GetMissedRatio()
{
	return totalMissed / ( totalPaired + totalFake );
}

//Return a root 2D histogram containing the smearing matrix
TH2F * SmearingMatrix::MakeRootHistogram( string Name, string Title, bool MakeNormalised )
{
	//Create the histogram - add some binwidths on for the underflow, overflow and junk bins
	int binNumber = indexCalculator->GetBinNumber() + 1;
	TH2F * outputHistogram = new TH2F( Name.c_str(), Title.c_str(), binNumber, 0.0, (double)binNumber, binNumber, 0.0, (double)binNumber );

	//Warn about bad method
	if ( !MakeNormalised && isFinalised)
	{
		cout << "WARNING: matrix has been normalised already, reversing this" << endl;
	}

	//Loop over all bins
	for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			int outputBin = outputHistogram->GetBin( causeIndex + 1, effectIndex + 1, 0 );

			//Deal with normalisation
			if ( MakeNormalised && !isFinalised )
			{
				outputHistogram->SetBinContent( outputBin, matrix[causeIndex][effectIndex] / normalisation[causeIndex] );
			}
			else if ( !MakeNormalised && isFinalised )
			{
				outputHistogram->SetBinContent( outputBin, matrix[causeIndex][effectIndex] * normalisation[causeIndex] );
			}
			else
			{
				outputHistogram->SetBinContent( outputBin, matrix[causeIndex][effectIndex] );
			}
		}
	}

	return outputHistogram;
}
