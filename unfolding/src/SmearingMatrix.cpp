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
SmearingMatrix::SmearingMatrix( Indices * InputIndices ) : indexCalculator(InputIndices), isFinalised(false)
{
	//Initialise the matrix and normalisation
	int binNumber = indexCalculator->GetBinNumber() + 1;

	vector<double> temporaryMatrixRow( binNumber, 0.0 );
	matrix = vector< vector<double> >( binNumber, temporaryMatrixRow );

	normalisation = vector<double>( binNumber, 0.0 );
}

//Populate the matrix with unnormalised input smearing matrices to combine
void SmearingMatrix::StoreUnnormalisedMatrix( TH2F * InputMatrix )
{
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
void SmearingMatrix::StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double Weight )
{
	//Look up the indices of the truth and reco
	int truthIndex = indexCalculator->GetIndex(Truth);
	int recoIndex = indexCalculator->GetIndex(Reco);

	//Increment values
	matrix[truthIndex][recoIndex] += Weight;
	normalisation[truthIndex] += Weight;
}
void SmearingMatrix::StoreUnreconstructedTruth( vector<double> Truth, double Weight )
{
	//Look up the index of the truth value
	int truthIndex = indexCalculator->GetIndex(Truth);
	int recoIndex = indexCalculator->GetBinNumber();

	//Increment values
	matrix[truthIndex][recoIndex] += Weight;
	normalisation[truthIndex] += Weight;
}
void SmearingMatrix::StoreReconstructedFake( vector<double> Reco, double Weight )
{
	//Look up the index of the reco value
	int truthIndex = indexCalculator->GetBinNumber();
	int recoIndex = indexCalculator->GetIndex(Reco);

	//Increment values
	matrix[truthIndex][recoIndex] += Weight;
	normalisation[truthIndex] += Weight;
}

//Do a bunch of extra calculations that aren't necessary if you just want the raw smearing matrix
void SmearingMatrix::Finalise()
{
	if ( !isFinalised )
	{
		//Normalise and calculate efficiencies
		int binNumber = indexCalculator->GetBinNumber();
		efficiencies = vector<double>( binNumber, 0.0 );
		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			//Calculate the fake rate for this effect index
			double fakeRate = matrix[binNumber][effectIndex] / binNumber;

			for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
			{
				//Add in the fakes, evenly distributed over all causes
				matrix[causeIndex][effectIndex] += fakeRate;

				//Fix nasty div0 errors caused by the treatment of fakes above
				if ( normalisation[causeIndex] == 0 && fakeRate == matrix[causeIndex][effectIndex] )
				{
					//Assume that the truth distribution has no entries in this bin, so it should purely have fakes
					matrix[causeIndex][effectIndex] = 1.0 / binNumber;
				}
				else
				{
					//Normalise
					matrix[causeIndex][effectIndex] /= normalisation[causeIndex];
				}

				//Check for divide-by-zero problems
				if ( isinf( matrix[causeIndex][effectIndex] ) )
				{
					cerr << "SMEARING: Divide by zero error (truth) in the normalisation" << endl;
					exit(1);
				}
				else if ( isnan( matrix[causeIndex][effectIndex] ) )
				{
					matrix[causeIndex][effectIndex] = 0.0;
				}

				efficiencies[causeIndex] += matrix[causeIndex][effectIndex];
			}
		}

		//Also normalise the fakes and unreconstructed events, just for nice presentation
		for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
		{
			//Normalise unreconstructed
			matrix[causeIndex][binNumber] /= normalisation[causeIndex];
		}
		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			//Normalise fakes
			matrix[binNumber][effectIndex] /= normalisation[binNumber];
		}

		//Initialise the tensor, apologies for the horrible syntax
		/*smearingCovariance = vector< vector< vector<double> > >( binNumber, vector< vector<double> >( binNumber, vector<double>( binNumber, 0.0 ) ) );

		//Loop over all r, s and u to populate the tensor
		for ( int indexR = 0; indexR < binNumber; indexR++ )
		{
		for ( int indexS = 0; indexS < binNumber; indexS++ )
		{
		for ( int indexU = 0; indexU < binNumber; indexU++ )
		{
		//Calculate the term
		double term = matrix[indexU][indexR];
		if ( indexR != indexS )
		{
		term *= -matrix[indexU][indexS];
		}
		else
		{
		term *= ( 1 - term );
		}

		term /= normalisation[indexU];

		//Check for divide-by-zero errors
		if ( isinf(term) )
		{
		cerr << "SMEARING: Divide by zero error (truth) in the error analysis" << endl;
		exit(1);
		}
		else if ( isnan(term) )
		{
		term = 0.0;
		}

		//Store the term
		smearingCovariance[indexR][indexS][indexU] = term;
		}
		}
		}*/

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

	return matrix[CauseIndex][EffectIndex];
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

//Return smearing covariance tensor element
double SmearingMatrix::GetCovarianceElement( int IndexR, int IndexS, int IndexU )
{
	cerr << "Don't use this. It's a waste of time" << endl;
	exit(1);

	//Check if smearing covariance is already calculated
	/*if ( !isFinalised )
	  {
	  Finalise();
	  }

	  return smearingCovariance[IndexR][IndexS][IndexU];*/
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
