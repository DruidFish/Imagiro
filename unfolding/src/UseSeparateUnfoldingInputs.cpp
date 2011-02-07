/**
  @class UseSeparateUnfoldingInputs

  The class that performs the unfolding by using the ROOT objects made by MakeSeparateUnfoldingInputs

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "UniformPrior.h"

#include "UseSeparateUnfoldingInputs.h"
#include "Distribution.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>

//Default constructor - useless
UseSeparateUnfoldingInputs::UseSeparateUnfoldingInputs()
{
}

//Constructor taking the required bin number,
//minimum and maximum of the output distribution as arguments
//NB: the unfolding scales roughly with bin number ^ 2, the
//error calculation scales roughly with bin number ^ 3.
UseSeparateUnfoldingInputs::UseSeparateUnfoldingInputs( int BinNumber, double Minimum, double Maximum, bool DebugMode ) : debug(DebugMode)
{
	indexCalculator = new Indices( vector<int>( 1, BinNumber ), vector<double>( 1, Minimum ), vector<double>( 1, Maximum ) );
	inputSmearing = new SmearingMatrix(indexCalculator);
}

//N-Dimensional version
UseSeparateUnfoldingInputs::UseSeparateUnfoldingInputs( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima, bool DebugMode ) : debug(DebugMode)
{
	indexCalculator = new Indices( BinNumbers, Minima, Maxima );
	inputSmearing = new SmearingMatrix(indexCalculator);
}

//Destructor
UseSeparateUnfoldingInputs::~UseSeparateUnfoldingInputs()
{
	delete indexCalculator;
	delete inputSmearing;
}

//If you want to make a smearing matrix from a bunch of
//other ones stored as (un-normalised) TH2Fs, use this
void UseSeparateUnfoldingInputs::StoreSmearingMatrix( TH2F * InputMatrix )
{
	inputSmearing->StoreUnnormalisedMatrix(InputMatrix);
}

//Store distribution(s) to use as a prior for the unfolding
void UseSeparateUnfoldingInputs::StorePriorDistribution( TH1F * InputPrior )
{
	inputPriorDistributions.push_back(InputPrior);
}

//Store data distribution(s) to combine and unfold
void UseSeparateUnfoldingInputs::StoreDataDistribution( TH1F * InputData )
{
	inputDataDistributions.push_back(InputData);
}

//Once all data is stored, run the unfolding
//You can specify when the iterations should end,
//with an upper limit on iteration number, or by
//comparing the results from the last two iterations.
//Iteration ends if the chi squared comparison value of
//the two results is lower than the threshold, or if
//the Kolmogorov-Smirnof comparison value is higher
void UseSeparateUnfoldingInputs::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	//Make the truth and data disributions
	Distribution * priorDistribution;
	Distribution * dataDistribution;
	if ( inputPriorDistributions.size() > 0 )
	{
		priorDistribution = new Distribution( inputPriorDistributions, indexCalculator );
	}
	else
	{
		cerr << "ERROR: No prior distribution provided" << endl;
		cerr << "This is the prior probability distribution in the Bayesian calculation" << endl;
		cerr << "If you want to use a uniform prior you could complain to the author (bwynne@cern.ch)" << endl;
		cerr << "Alternatively there is a UniformPrior class that you could use - edit the code" << endl;
		exit(1);
	}
	if ( inputDataDistributions.size() > 0 )
	{
		dataDistribution = new Distribution( inputDataDistributions, indexCalculator );
	}
	else
	{
		cerr << "ERROR: No data distribution provided" << endl;
		cerr << "Either you're trying to make magic data out of thin air or you've forgotten something..." << endl;
		exit(1);
	}

	//Make a histogram of the truth distribution
	TH1F * priorHistogram = priorDistribution->MakeRootHistogram( "truth", "Monte Carlo truth distribution" );
	allResults.push_back(priorHistogram);

	//Normalise the smearing matrix
	inputSmearing->Finalise();

	//Debug output
	TFile * debugFile;
	if (debug)
	{
		debugFile = new TFile( "unfoldingDebugOutput.root", "RECREATE" );
		priorHistogram->Write();
		dataDistribution->MakeRootHistogram( "data", "Uncorrected data distribution" )->Write();
		inputSmearing->MakeRootHistogram( "smearing", "Smearing matrix" )->Write();
	}

	//Iterate, making new distribution from data, old distribution and smearing matrix
	Distribution * adjustedDistribution;
	TH1F * adjustedHistogram;
	char name[20];
	for ( int iteration = 0; iteration < MostIterations; iteration++ )
	{
		//Iterate
		adjustedDistribution = new Distribution( dataDistribution, inputSmearing, priorDistribution, indexCalculator );
		if (WithSmoothing)
		{
			adjustedDistribution->Smooth();
		}

		//Make a root histogram
		sprintf( name, "iteration%d", iteration );
		adjustedHistogram = adjustedDistribution->MakeRootHistogram( name, name );
		allResults.push_back(adjustedHistogram);

		//Compare with previous distribution
		double chi2 = adjustedHistogram->Chi2Test( priorHistogram, "UUCHI2" );
		double kolmogorov = adjustedHistogram->KolmogorovTest( priorHistogram, "" );

		//Replace previous iteration
		priorDistribution = adjustedDistribution;
		priorHistogram = adjustedHistogram;

		//Debug output
		if (debug)
		{
			cout << "Chi squared = " << chi2 << " and K-S probability = " << kolmogorov << " for iteration " << iteration << endl;
			adjustedHistogram->Write();
		}

		//Check for termination conditions
		if ( chi2 < ChiSquaredThreshold || kolmogorov > KolmogorovThreshold )
		{
			cout << "Converged" << endl;
			lastResult = adjustedDistribution;
			break;
		}
		else if ( iteration == MostIterations - 1 )
		{
			cerr << "No convergence after " << MostIterations << " iterations" << endl;
			lastResult = adjustedDistribution;
			break;
		}
	}

	//Close the debug file
	if (debug)
	{
		debugFile->Close();
	}
}

//Retrieve a TH1F* containing the unfolded data
//distribution, with or without errors
//NB: the error calculation is only performed
//when you run the method with errors for the first time
TH1F * UseSeparateUnfoldingInputs::UnfoldedDistribution( bool WithErrors )
{
	if (WithErrors)
	{
		return lastResult->MakeRootHistogram( "unfolded", "Unfolded distribution with errors", true );
	}
	else
	{
		return lastResult->MakeRootHistogram( "unfolded", "Unfolded distribution", false );
	}
}

//Retrieve the unfolded distribution from each iteration
//as a vector of TH1F*.
//The 0th entry will be the prior (MC truth) distirbution,
//and the last entry will be the same as that returned
//by UnfoldedDistribution()
vector< TH1F* > UseSeparateUnfoldingInputs::AllIterationResults()
{
	return allResults;
}
