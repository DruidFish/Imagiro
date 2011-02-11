/**
  @class IterativeUnfolding

  The class that performs unfolding on a data distribution in the simple case when all events are available

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "IterativeUnfolding.h"
#include "Distribution.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>

const int MAX_ITERATIONS_FOR_CROSS_CHECK = 10;

//Default constructor - useless
IterativeUnfolding::IterativeUnfolding()
{
}

//Constructor taking the required bin number,
//minimum and maximum of the output distribution as arguments
//NB: the unfolding scales roughly with bin number ^ 2, the
//error calculation scales roughly with bin number ^ 3.
IterativeUnfolding::IterativeUnfolding( int BinNumber, double Minimum, double Maximum, string Name, int UniqueID, bool DebugMode ) : debug(DebugMode), name(Name), uniqueID(UniqueID)
{
	//Make new vectors just with the one entry
	indexCalculator = new Indices( vector<int>( 1, BinNumber ), vector<double>( 1, Minimum ), vector<double>( 1, Maximum ) );
	inputSmearing = new SmearingMatrix(indexCalculator);
	priorDistribution = new Distribution(indexCalculator);
	dataDistribution = new Distribution(indexCalculator);
	simulatedDistribution = new Distribution(indexCalculator);
	sumOfDataWeightSquares = vector<double>( indexCalculator->GetBinNumber(), 0.0 );
}

//N-Dimensional version
IterativeUnfolding::IterativeUnfolding( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima, string Name, int UniqueID, bool DebugMode ) : debug(DebugMode), name(Name), uniqueID(UniqueID)
{
	indexCalculator = new Indices( BinNumbers, Minima, Maxima );
	inputSmearing = new SmearingMatrix(indexCalculator);
	priorDistribution = new Distribution(indexCalculator);
	dataDistribution = new Distribution(indexCalculator);
	simulatedDistribution = new Distribution(indexCalculator);
	sumOfDataWeightSquares = vector<double>( indexCalculator->GetBinNumber(), 0.0 );
}

//Destructor
IterativeUnfolding::~IterativeUnfolding()
{
	delete indexCalculator;
	delete inputSmearing;
	delete priorDistribution;
	delete dataDistribution;
	delete simulatedDistribution;
	allResults.clear();
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void IterativeUnfolding::StoreTruthRecoPair( double Truth, double Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if (UseInPrior)
	{
		priorDistribution->StoreEvent( vector<double>( 1, Truth ), TruthWeight );
		simulatedDistribution->StoreEvent( vector<double>( 1, Reco ), RecoWeight );
	}
	inputSmearing->StoreTruthRecoPair( vector<double>( 1, Truth ), vector<double>( 1, Reco ), TruthWeight * RecoWeight );
}

//N-Dimensional version
void IterativeUnfolding::StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if (UseInPrior)
	{
		priorDistribution->StoreEvent( Truth, TruthWeight );
		simulatedDistribution->StoreEvent( Reco, RecoWeight );
	}
	inputSmearing->StoreTruthRecoPair( Truth, Reco, TruthWeight * RecoWeight );
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void IterativeUnfolding::StoreUnreconstructedTruth( double Truth, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		priorDistribution->StoreEvent( vector<double>( 1, Truth ), Weight );
	}
	inputSmearing->StoreUnreconstructedTruth( vector<double>( 1, Truth ), Weight );
}

//N-Dimensional version
void IterativeUnfolding::StoreUnreconstructedTruth( vector<double> Truth, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		priorDistribution->StoreEvent( Truth, Weight );
	}
	inputSmearing->StoreUnreconstructedTruth( Truth, Weight );
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void IterativeUnfolding::StoreReconstructedFake( double Reco, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		simulatedDistribution->StoreEvent( vector<double>( 1, Reco ), Weight );
	}
	inputSmearing->StoreReconstructedFake( vector<double>( 1, Reco ), Weight );
}

//N-Dimensional version
void IterativeUnfolding::StoreReconstructedFake( vector<double> Reco, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		simulatedDistribution->StoreEvent( Reco, Weight );
	}
	inputSmearing->StoreReconstructedFake( Reco, Weight );
}

//Store a value from the uncorrected data distribution
void IterativeUnfolding::StoreDataValue( double Data, double Weight )
{
	vector<double> dataVector( 1, Data );
	dataDistribution->StoreEvent( dataVector, Weight );
	sumOfDataWeightSquares[ indexCalculator->GetIndex( dataVector ) ] += ( Weight * Weight );
}

//N-Dimensional version
void IterativeUnfolding::StoreDataValue( vector<double> Data, double Weight )
{
	dataDistribution->StoreEvent( Data, Weight );
	sumOfDataWeightSquares[ indexCalculator->GetIndex( Data ) ] += ( Weight * Weight );
}

//Once all data is stored, run the unfolding
//You can specify when the iterations should end,
//with an upper limit on iteration number, or by
//comparing the results from the last two iterations.
//Iteration ends if the chi squared comparison value of
//the two results is lower than the threshold, or if
//the Kolmogorov-Smirnof comparison value is higher
void IterativeUnfolding::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	//Make a histogram of the truth distribution
	char plotName[ name.size() + 20 ];
	sprintf( plotName, "%sPrior%d", name.c_str(), uniqueID );
	TH1F * priorHistogram = priorDistribution->MakeRootHistogram( plotName, "Monte Carlo prior distribution" );
	allResults.push_back(priorHistogram);

	//Make another histogram to save for use later
	sprintf( plotName, "%sTruth%d", name.c_str(), uniqueID );
	truthHistogram = priorDistribution->MakeRootHistogram( plotName, "Monte Carlo truth distribution" );

	//Finalise the smearing matrix
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
	char iterationName[ name.size() + 20 ];
	for ( int iteration = 0; iteration < MostIterations; iteration++ )
	{
		//Smooth the prior distribution
		if ( WithSmoothing )
		{
			priorDistribution->Smooth();
		}

		//Iterate
		adjustedDistribution = new Distribution( dataDistribution, inputSmearing, priorDistribution, indexCalculator );

		//Make a root histogram
		sprintf( iterationName, "%s%dIteration%d", name.c_str(), uniqueID, iteration );
		adjustedHistogram = adjustedDistribution->MakeRootHistogram( iterationName, iterationName );
		allResults.push_back(adjustedHistogram);

		//Compare with previous distribution
		double chi2 = adjustedHistogram->Chi2Test( priorHistogram, "UUCHI2" );
		double kolmogorov = adjustedHistogram->KolmogorovTest( priorHistogram, "" );

		//Reset for next iteration
		delete priorDistribution;
		priorDistribution = adjustedDistribution;
		priorHistogram = adjustedHistogram;

		//Debug output
		if (debug)
		{
			cout << "Chi squared = " << chi2 << " and K-S probability = " << kolmogorov << " for iteration " << iteration << endl;
			adjustedHistogram->Write();
		}

		//Check for termination conditions
		if ( chi2 < ChiSquaredThreshold )
		{
			cout << "Converged: chi2 " << chi2 << " < " << ChiSquaredThreshold << endl;
			lastResult = adjustedDistribution;
			break;
		}
		if ( kolmogorov > KolmogorovThreshold )
		{
			cout << "Converged: KS " << kolmogorov << " > " << KolmogorovThreshold << endl;
			lastResult = adjustedDistribution;
			break;
		}
		if ( iteration == MostIterations - 1 )
		{
			cerr << "Maximum number of iterations (" << MostIterations << ") reached" << endl;
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

//Perform a closure test
//Unfold the MC reco distribution with the corresponding truth information as a prior
//It should give the truth information back exactly...
void IterativeUnfolding::ClosureTest()
{
	//Make a histogram for the Root comparisons
	char plotName[ name.size() + 20 ];
	sprintf( plotName, "%sClosurePrior%d", name.c_str(), uniqueID );
	TH1F * priorHistogram = priorDistribution->MakeRootHistogram( plotName, "Closure test prior distribution" );

	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Unfold once only
	Distribution * adjustedDistribution = new Distribution( simulatedDistribution, inputSmearing, priorDistribution, indexCalculator );
	sprintf( plotName, "%sClosureResult%d", name.c_str(), uniqueID );
	TH1F * adjustedHistogram = adjustedDistribution->MakeRootHistogram( plotName, plotName );

	//Compare with truth distribution
	double chi2 = adjustedHistogram->Chi2Test( priorHistogram, "UUCHI2" );
	double kolmogorov = adjustedHistogram->KolmogorovTest( priorHistogram, "" );
	cout << endl << "Closure test - comparing unfolded MC reco with MC truth: " << endl;
	cout << "Chi squared = " << chi2 << " and K-S probability = " << kolmogorov << endl;
}

//Perform an unfolding cross-check
//Use MC truth A as a prior to unfold MC reco B
//Iterations cease when result is sufficiently close to MC truth B (passed as argument)
//Returns the number of iterations required. Convergence criteria as output arguments
int IterativeUnfolding::MonteCarloCrossCheck( TH1F * ReferencePlot, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	//Make a histogram of the truth distribution
	char plotName[ name.size() + 20 ];
	sprintf( plotName, "%sRawPrior%d", name.c_str(), uniqueID );
	TH1F * priorHistogram = priorDistribution->MakeRootHistogram( plotName, "Monte Carlo raw prior distribution" );

	//Make another histogram to save for use later
	sprintf( plotName, "%sTruth%d", name.c_str(), uniqueID );
	truthHistogram = priorDistribution->MakeRootHistogram( plotName, "Monte Carlo truth distribution" );

	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Compare the uncorrected data to the truth
	TH1F * dataHistogram = dataDistribution->MakeRootHistogram( plotName, "Uncorrected data distribution for cross check" );
	double lastChiSquared = dataHistogram->Chi2Test( ReferencePlot, "UUCHI2" );
	double lastKolmogorov = dataHistogram->KolmogorovTest( ReferencePlot, "" );
	double chiSquaredResult = lastChiSquared;
	double kolmogorovResult = lastKolmogorov;
	delete dataHistogram;

	//Iterate, making new distribution from data, old distribution and smearing matrix
	Distribution * adjustedDistribution;
	TH1F * adjustedHistogram;
	char iterationName[ name.size() + 20 ];
	for ( int iteration = 0; iteration < MAX_ITERATIONS_FOR_CROSS_CHECK; iteration++ )
	{
		//Smooth the prior distribution
		if ( WithSmoothing )
		{
			priorDistribution->Smooth();
		}

		//Iterate
		adjustedDistribution = new Distribution( dataDistribution, inputSmearing, priorDistribution, indexCalculator );

		//Make a root histogram for comparison
		sprintf( iterationName, "%s%dIteration%d", name.c_str(), uniqueID, iteration );
		adjustedHistogram = adjustedDistribution->MakeRootHistogram( iterationName, iterationName );

		//Compare with reference distribution
		double referenceChi2 = adjustedHistogram->Chi2Test( ReferencePlot, "UUCHI2" );
		double referenceKolmogorov = adjustedHistogram->KolmogorovTest( ReferencePlot, "" );

		//Compare with last iteration
		double chi2 = adjustedHistogram->Chi2Test( priorHistogram, "UUCHI2" );
		double kolmogorov = adjustedHistogram->KolmogorovTest( priorHistogram, "" );

		//Reset for next iteration
		delete priorDistribution;
		delete priorHistogram;
		priorDistribution = adjustedDistribution;
		priorHistogram = adjustedHistogram;

		//Check to see if things have got worse
		if ( referenceChi2 > lastChiSquared || referenceKolmogorov < lastKolmogorov || iteration == MAX_ITERATIONS_FOR_CROSS_CHECK - 1 )
		{
			//Return the criteria
			ChiSquaredThreshold = chiSquaredResult;
			KolmogorovThreshold = kolmogorovResult;
			return iteration;
		}
		else
		{
			//Update the last values
			lastChiSquared = referenceChi2;
			lastKolmogorov = referenceKolmogorov;
			chiSquaredResult = chi2;
			kolmogorovResult = kolmogorov;
		}
	}
}

//Retrieve a TH1F* containing the unfolded data
//distribution, with or without errors
//NB: the error calculation is only performed
//when you run the method with errors for the first time
TH1F * IterativeUnfolding::UnfoldedDistribution( string Name, string Title, bool WithErrors )
{
	if (WithErrors)
	{
		Title += " with errors";
		return lastResult->MakeRootHistogram( Name, Title, true );
	}
	else
	{
		return lastResult->MakeRootHistogram( Name, Title, false );
	}
}

//Retrieve the unfolded distribution from each iteration
//as a vector of TH1F*.
//The 0th entry will be the prior (MC truth) distirbution,
//and the last entry will be the same as that returned
//by UnfoldedDistribution()
vector< TH1F* > IterativeUnfolding::AllIterationResults()
{
	return allResults;
}

//Retrieve the smearing matrix used
TH2F * IterativeUnfolding::GetSmearingMatrix( string Name, string Title )
{
	return inputSmearing->MakeRootHistogram( Name, Title );
}

//Retrieve the truth distribution
TH1F * IterativeUnfolding::GetTruthDistribution( string Name, string Title )
{
	//Check that the histogram has been made (since we need it before the unfolding starts)
	if ( !truthHistogram )
	{
		truthHistogram = priorDistribution->MakeRootHistogram( Name, Title );
	}
	else
	{
		truthHistogram->SetName( Name.c_str() );
		truthHistogram->SetTitle( Title.c_str() );
	}
	return truthHistogram;
}

//Retrieve the uncorrected data distribution
TH1F * IterativeUnfolding::GetUncorrectedDataDistribution( string Name, string Title )
{
	return dataDistribution->MakeRootHistogram( Name, Title );
}

//Handy for error calculation
vector<double> IterativeUnfolding::SumOfDataWeightSquares()
{
	return sumOfDataWeightSquares;
}
