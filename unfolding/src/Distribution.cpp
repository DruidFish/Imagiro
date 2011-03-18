/**
  @class Distribution

  A 1D data histogram / probability distribution. Contains the folding and unfolding methods

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "Distribution.h"
#include "UnfoldingMatrix.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

//Default constructor
Distribution::Distribution()
{
}

//Just initialise an empty distribution
Distribution::Distribution( Indices * InputIndices ) : indexCalculator(InputIndices)
{
	//Initialise the bins (include a bad bin)
	binValues = vector<double>( InputIndices->GetBinNumber() + 1, 0.0 );
	integral = 0.0;
}

//Initialise by summing a bunch of TH1Fs
Distribution::Distribution( vector< TH1F* > InputDistributions, Indices * InputIndices ) : integral(0.0), indexCalculator(InputIndices)
{
	//Get the number of bins in the distribution (include a bad bin)
	int binNumber = InputIndices->GetBinNumber() + 1;

	//Check all inputs have the right number of bins
	for ( int inputIndex = 0; inputIndex < InputDistributions.size(); inputIndex++ )
	{
		int inputBinNumber = InputDistributions[inputIndex]->GetNbinsX();
		if ( inputBinNumber != binNumber - 2 )
		{
			cerr << "ERROR: input data distribution does not have the correct number of bins: " << inputBinNumber << " vs " << binNumber << endl;
			exit(1);
		}
	}

	//Initialise and populate the distribution
	binValues = vector<double>( binNumber, 0.0 );
	for ( int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		for ( int inputIndex = 0; inputIndex < InputDistributions.size(); inputIndex++ )
		{
			binValues[binIndex] += InputDistributions[inputIndex]->GetBinContent(binIndex);
			integral += InputDistributions[inputIndex]->GetBinContent(binIndex);
		}
	}
}

//Calculate a corrected distribtion
Distribution::Distribution( Distribution * DataDistribution, SmearingMatrix * Smearing, Distribution * PriorDistribution ) : integral(0.0),
	indexCalculator( DataDistribution->indexCalculator )
{
	//Get the number of bins in the distribution (include a bad bin)
	int binNumber = indexCalculator->GetBinNumber() + 1;

	//Make a new, empty distribution
	binValues = vector< double >( binNumber, 0.0 );

	//Create the Bayes theorem inverse of the smearing matrix
	UnfoldingMatrix * bayes = new UnfoldingMatrix( Smearing, PriorDistribution, indexCalculator );

	//Populate the distribution
	for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
	{
		//Retrieve the uncorrected value
		double effectNumber = DataDistribution->GetBinNumber( effectIndex );

		//Apply the unfolding matrix
		for  ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
		{
			double newValue = bayes->GetElement( causeIndex, effectIndex ) * effectNumber;
			binValues[ causeIndex ] += newValue;
			integral += newValue;
		}
	}

	//Free up some memory
	delete bayes;
}

//Make this distribution by smearing another
Distribution::Distribution( Distribution * InputDistribution, SmearingMatrix * Smearing ) : integral(0.0), indexCalculator( InputDistribution->indexCalculator )
{
	//Get the number of bins in the distribution (include a bad bin)
	int binNumber = indexCalculator->GetBinNumber() + 1;

	//Make a new, empty distribution
	binValues = vector< double >( binNumber, 0.0 );

	//Loop over each bin in this distribution
	for ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		//Retrieve the correct value
		double causeNumber = InputDistribution->GetBinNumber( causeIndex );

		//Apply the smearing matrix
		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			double newValue = Smearing->GetElement( causeIndex, effectIndex ) * causeNumber;
			binValues[ effectIndex ] += newValue;
			integral += newValue;
		}
	}
}

//Destructor
Distribution::~Distribution()
{
}

//Store an event
void Distribution::StoreEvent( vector<double> Value, double Weight )
{
	int binIndex = indexCalculator->GetIndex( Value );
	binValues[ binIndex ] += Weight;
	integral += Weight;
}

void Distribution::StoreBadEvent( double Weight )
{
	binValues[ binValues.size() - 1 ] += Weight;
	integral += Weight;
}

void Distribution::SetBadBin( double Ratio )
{
	//Check for existing bad bin values
	if ( binValues[ binValues.size() - 1 ] != 0.0 )
	{
		cerr << "WARNING: Overwriting bad bin value suggests that unnecessary extrapolation is being made" << endl;
		integral -= binValues[ binValues.size() - 1 ];
	}

	//Set the new bin value as the total volume of the distribution scaled by the given ratio
	binValues[ binValues.size() - 1 ] = integral * Ratio;
	integral += binValues[ binValues.size() - 1 ];
}

//Return the contents of the bin with the given index
double Distribution::GetBinNumber( int InputIndex )
{
	if ( InputIndex >= 0 && InputIndex <= indexCalculator->GetBinNumber() )
	{
		return binValues[InputIndex];
	}
	else
	{
		cerr << "FATAL: Invalid bin index " << InputIndex << endl;
		exit(1);
	}
}

//Return the normalised contents of the bin with the given index
double Distribution::GetBinProbability( int InputIndex )
{
	if ( InputIndex >= 0 && InputIndex <= indexCalculator->GetBinNumber() )
	{
		return (double)( binValues[InputIndex] ) / (double)integral;
	}
	else
	{
		cerr << "FATAL: Invalid bin index " << InputIndex << endl;
		exit(1);
	}
}

//Make a root histogram from the distribution
TH1F * Distribution::MakeRootHistogram( string Name, string Title, bool MakeNormalised, bool WithBadBin )
{
	//Get the correct number of bins
	int binNumber = indexCalculator->GetBinNumber();
	if (WithBadBin)
	{
		binNumber += 1;
	}

	//Make a new root histogram
	TH1F * rootHistogram = new TH1F( Name.c_str(), Title.c_str(), binNumber - 2, indexCalculator->GetMinima()[0], indexCalculator->GetMaxima()[0] );

	//Populate the bins one-by-one - note that in the TH1F, bin 0 is the underflow
	for ( int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		//Normalise the distribution (or not)
		if (MakeNormalised)
		{
			rootHistogram->SetBinContent( binIndex, (double)( binValues[binIndex] ) / (double)integral );
		}
		else
		{
			rootHistogram->SetBinContent( binIndex, binValues[binIndex] );
		}
	}

	return rootHistogram;
}

//Smooth the distribution using moving average
void Distribution::Smooth( int SideBinNumber )
{
	int binNumber = indexCalculator->GetBinNumber();
	vector<double> newBinValues;
	vector<int> separateBinIndices, separateSumIndices;
	double newIntegral = 0.0;

	//Operate on all bins
	for ( int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		//Don't smooth the over/underflow bins of the first dimension
		separateBinIndices = indexCalculator->GetNDimensionalIndex( binIndex );
		if ( separateBinIndices[0] == 0 || separateBinIndices[0] == indexCalculator->GetBinNumber(0) - 1 )
		{
			//Just copy the original value
			newBinValues.push_back( binValues[binIndex] );
			newIntegral += binValues[binIndex];
		}
		else
		{
			double binTotal = 0.0;

			//Take the average of the local bins
			for ( int localIndex = -SideBinNumber; localIndex <= SideBinNumber; localIndex++ )
			{
				int sumIndex = binIndex + localIndex;

				//Can this bin be used in averaging?
				bool useInAveraging = true;
				separateSumIndices = indexCalculator->GetNDimensionalIndex( sumIndex );
				if ( separateSumIndices[0] == 0 || separateSumIndices[0] == indexCalculator->GetBinNumber(0) - 1 )
				{
					useInAveraging = false;
				}
				else
				{
					for ( int dimensionIndex = 1; dimensionIndex < separateBinIndices.size(); dimensionIndex++ )
					{
						if ( separateSumIndices[dimensionIndex] != separateBinIndices[dimensionIndex] )
						{
							//Bins aren't in the same part of the distribution
							useInAveraging = false;
							break;
						}
					}
				}

				if ( sumIndex >= 0 && sumIndex < binNumber && useInAveraging )
				{
					//Use this bin in the averaging
					binTotal += binValues[sumIndex];
				}
				else
				{
					//Derive an edge-preserving value instead
					binTotal += ( 2.0 * binValues[binIndex] ) - binValues[binIndex - localIndex];
				}
			}

			//Store the average as a new distribution
			binTotal /= ( ( 2.0 * (double)SideBinNumber ) + 1 );
			newBinValues.push_back(binTotal);
			newIntegral += binTotal;
		}
	}

	binValues = newBinValues;
	integral = newIntegral;
}
