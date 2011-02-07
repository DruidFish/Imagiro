/**
  @class Distribution

  A 1D data histogram / probability distribution. Contains the unfolding method

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "Distribution.h"
#include "UnfoldingMatrix.h"
//#include "CovarianceMatrix.h"
//#include "JustVariances.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

//Default constructor
Distribution::Distribution()
{
}

//Just initialise directly from a list a data points
/*Distribution::Distribution( vector<double> InputData, Indices * InputIndices ) : indexCalculator(InputIndices)
  {
//Initialise the bins
int binNumber = InputIndices->GetBinNumber();
binValues = vector<double>( binNumber, 0.0 );
integral = InputData.size();

//Populate the distribution
for ( int eventIndex = 0; eventIndex < InputData.size(); eventIndex++ )
{
//Look up the index of the data point
int binIndex = InputIndices->GetIndex( InputData[eventIndex] );

//Add to the appropriate values
if ( binIndex >= 0 )
{
binValues[binIndex]++;
}
}
}*/

//Just initialise an empty distribution
Distribution::Distribution( Indices * InputIndices ) : indexCalculator(InputIndices)
{
	//Initialise the bins
	binValues = vector<double>( InputIndices->GetBinNumber(), 0.0 );
	integral = 0.0;
}

//Initialise by summing a bunch of TH1Fs
Distribution::Distribution( vector< TH1F* > InputDistributions, Indices * InputIndices ) : indexCalculator(InputIndices)
{
	//Check all inputs have the right number of bins
	int binNumber = InputIndices->GetBinNumber();
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
	integral = 0.0;
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
Distribution::Distribution( Distribution * DataDistribution, SmearingMatrix * Smearing, Distribution * PriorDistribution, Indices * InputIndices ) : integral(0.0),
	indexCalculator(InputIndices)
{
	//Create the Bayes theorem inverse of the smearing matrix
	UnfoldingMatrix * bayes = new UnfoldingMatrix( Smearing, PriorDistribution, InputIndices );

	//Populate the distribution
	int binNumber = InputIndices->GetBinNumber();
	for  ( int causeIndex = 0; causeIndex < binNumber; causeIndex++ )
	{
		double binValue = 0.0;

		for ( int effectIndex = 0; effectIndex < binNumber; effectIndex++ )
		{
			binValue += ( DataDistribution->GetBinNumber(effectIndex) * bayes->GetElement( causeIndex, effectIndex ) );
		}

		binValues.push_back(binValue);
		integral += binValue;
	}

	//Make the covariance matrix
	//errorCalculator = new JustVariances( bayes, Smearing, DataDistribution, integral, InputIndices );

	//Free up some memory
	delete bayes;
}

//Destructor
Distribution::~Distribution()
{
}

//Store an event
void Distribution::StoreEvent( vector<double> Value, double Weight )
{
	int binIndex = indexCalculator->GetIndex(Value);
	binValues[binIndex] += Weight;
	integral += Weight;
}

//Return the contents of the bin with the given index
double Distribution::GetBinNumber( int InputIndex )
{
	if ( InputIndex >= 0 && InputIndex < indexCalculator->GetBinNumber() )
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
	if ( InputIndex >= 0 && InputIndex < indexCalculator->GetBinNumber() )
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
TH1F * Distribution::MakeRootHistogram( string Name, string Title, bool WithErrors, bool MakeNormalised )
{
	int binNumber = indexCalculator->GetBinNumber();
	TH1F * distributionHistogram = new TH1F( Name.c_str(), Title.c_str(), binNumber - 2, indexCalculator->GetMinima()[0], indexCalculator->GetMaxima()[0] );

	//Populate the bins one-by-one - note that in the TH1F, bin 0 is the underflow
	for ( int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		//Normalise the distribution (or not)
		if (MakeNormalised)
		{
			distributionHistogram->SetBinContent( binIndex, (double)( binValues[binIndex] ) / (double)integral );
		}
		else
		{
			distributionHistogram->SetBinContent( binIndex, binValues[binIndex] );
		}

		//Add the errors if required
		//if ( WithErrors && errorCalculator )
		//{
		//	distributionHistogram->SetBinError( binIndex, errorCalculator->GetStandardDeviation(binIndex) );
		//}
	}

	return distributionHistogram;
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
