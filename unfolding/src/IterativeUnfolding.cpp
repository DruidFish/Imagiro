/**
  @class IterativeUnfolding

  The class that manages unfolding a distribution, and tests the quality of the unfolding

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "IterativeUnfolding.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <sstream>

const unsigned int MAX_ITERATIONS_FOR_CROSS_CHECK = 10;

//Default constructor - useless
IterativeUnfolding::IterativeUnfolding()
{
}

//Constructor taking the required bin number,
//minimum and maximum of the output distribution as arguments
//NB: the unfolding scales roughly with bin number ^ 2, the
//error calculation scales roughly with bin number ^ 3.
IterativeUnfolding::IterativeUnfolding( unsigned int BinNumber, double Minimum, double Maximum, string Name, unsigned int UniqueID, bool DebugMode )
{
	debug = DebugMode;
	name = Name;
	uniqueID = UniqueID;
	totalPaired = 0.0;
	totalFake = 0.0;
	totalMissed = 0.0;

	indexCalculator = new Indices( vector< unsigned int >( 1, BinNumber ), vector< double >( 1, Minimum ), vector< double >( 1, Maximum ) );

	inputSmearing = new SmearingMatrix( indexCalculator );

	dataDistribution = new Distribution( indexCalculator );
	unfoldedDistribution = new Distribution( indexCalculator );
	truthDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = new Distribution( indexCalculator );

	sumOfDataWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );

	distributionComparison = new Comparison( Name, UniqueID );
}

//N-Dimensional version
IterativeUnfolding::IterativeUnfolding( vector< unsigned int > BinNumbers, vector< double > Minima, vector< double > Maxima, string Name, unsigned int UniqueID, bool DebugMode )
{
	debug = DebugMode;
	name = Name;
	uniqueID = UniqueID;
	totalPaired = 0.0;
	totalFake = 0.0;
	totalMissed = 0.0;

	indexCalculator = new Indices( BinNumbers, Minima, Maxima );

	inputSmearing = new SmearingMatrix( indexCalculator );

	dataDistribution = new Distribution( indexCalculator );
	unfoldedDistribution = new Distribution( indexCalculator );
	truthDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = new Distribution( indexCalculator );

	sumOfDataWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );

	distributionComparison = new Comparison( Name, UniqueID );
}

//Destructor
IterativeUnfolding::~IterativeUnfolding()
{
	delete indexCalculator;

	delete inputSmearing;

	delete dataDistribution;
	delete unfoldedDistribution;
	delete truthDistribution;
	delete reconstructedDistribution;

	delete distributionComparison;
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void IterativeUnfolding::StoreTruthRecoPair( double Truth, double Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreEvent( vector< double >( 1, Truth ), TruthWeight );
		reconstructedDistribution->StoreEvent( vector< double >( 1, Reco ), RecoWeight );

	}

	totalPaired += TruthWeight;
	inputSmearing->StoreTruthRecoPair( vector< double >( 1, Truth ), vector< double >( 1, Reco ), TruthWeight, RecoWeight );
}

//N-Dimensional version
void IterativeUnfolding::StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreEvent( Truth, TruthWeight );
		reconstructedDistribution->StoreEvent( Reco, RecoWeight );

	}

	totalPaired += TruthWeight;
	inputSmearing->StoreTruthRecoPair( Truth, Reco, TruthWeight, RecoWeight );
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void IterativeUnfolding::StoreUnreconstructedTruth( double Truth, double Weight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreEvent( vector< double >( 1, Truth ), Weight );
		reconstructedDistribution->StoreBadEvent( Weight );
	}

	totalMissed += Weight;
	inputSmearing->StoreUnreconstructedTruth( vector< double >( 1, Truth ), Weight );
}

//N-Dimensional version
void IterativeUnfolding::StoreUnreconstructedTruth( vector< double > Truth, double Weight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreEvent( Truth, Weight );
		reconstructedDistribution->StoreBadEvent( Weight );
	}

	totalMissed += Weight;
	inputSmearing->StoreUnreconstructedTruth( Truth, Weight );
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void IterativeUnfolding::StoreReconstructedFake( double Reco, double Weight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreBadEvent( Weight );
		reconstructedDistribution->StoreEvent( vector< double >( 1, Reco ), Weight );
	}

	totalFake += Weight;
	inputSmearing->StoreReconstructedFake( vector< double >( 1, Reco ), Weight );
}

//N-Dimensional version
void IterativeUnfolding::StoreReconstructedFake( vector< double > Reco, double Weight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreBadEvent( Weight );
		reconstructedDistribution->StoreEvent( Reco, Weight );
	}

	totalFake += Weight;
	inputSmearing->StoreReconstructedFake( Reco, Weight );
}

//Store a value from the uncorrected data distribution
void IterativeUnfolding::StoreDataValue( double Data, double Weight )
{
	vector< double > dataVector( 1, Data );
	dataDistribution->StoreEvent( dataVector, Weight );
	sumOfDataWeightSquares[ indexCalculator->GetIndex( dataVector ) ] += ( Weight * Weight );
}

//N-Dimensional version
void IterativeUnfolding::StoreDataValue( vector< double > Data, double Weight )
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
void IterativeUnfolding::Unfold( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, unsigned int ErrorMode, bool WithSmoothing )
{
	//Extrapolate the number of missed events in the data
	dataDistribution->SetBadBin( totalMissed / ( totalPaired + totalFake ) );

	//Use the truth distribution as the prior
	Distribution * priorDistribution = truthDistribution;

	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Debug output
	TFile * debugFile;
	if ( debug )
	{
		debugFile = new TFile( "unfoldingDebugOutput.root", "RECREATE" );
		priorDistribution->MakeRootHistogram( "prior", "MC truth distribution as prior" )->Write();
		dataDistribution->MakeRootHistogram( "data", "Uncorrected data distribution" )->Write();
		inputSmearing->MakeRootHistogram( "smearing", "Smearing matrix" )->Write();
	}

	//Iterate, making new distribution from data, old distribution and smearing matrix
	UnfoldingMatrix * lastUnfoldingMatrix;
	for ( unsigned int iteration = 0; iteration < MostIterations; iteration++ )
	{
		//Smooth the prior distribution, if asked. Don't smooth the truth
		if ( WithSmoothing && iteration != 0 )
		{
			priorDistribution->Smooth();
		}

		//Make a new unfolding matrix
		if ( iteration != 0 )
		{
			delete lastUnfoldingMatrix;
		}
		lastUnfoldingMatrix = new UnfoldingMatrix( inputSmearing, priorDistribution );

		//Unfold
		unfoldedDistribution = new Distribution( dataDistribution, lastUnfoldingMatrix );

		//Compare with previous distribution
		double chi2, kolmogorov;
		distributionComparison->CompareDistributions( unfoldedDistribution, priorDistribution, chi2, kolmogorov, false );

		//Reset for next iteration
		if ( iteration != 0 )
		{
			//Don't delete the MC truth
			delete priorDistribution;
		}
		priorDistribution = unfoldedDistribution;

		//Debug output
		if ( debug )
		{
			cout << "Chi squared = " << chi2 << " and K-S probability = " << kolmogorov << " for iteration " << iteration << endl;
			stringstream iterationName;
			iterationName << "iteration" << iteration;
			unfoldedDistribution->MakeRootHistogram( iterationName.str(), iterationName.str() )->Write();
		}

		//Check for termination conditions
		if ( chi2 < ChiSquaredThreshold )
		{
			cout << "Converged: chi2 " << chi2 << " < " << ChiSquaredThreshold << endl;
			break;
		}
		if ( kolmogorov > KolmogorovThreshold )
		{
			cout << "Converged: KS " << kolmogorov << " > " << KolmogorovThreshold << endl;
			break;
		}
		if ( iteration == MostIterations - 1 )
		{
			cout << "Maximum number of iterations (" << MostIterations << ") reached" << endl;
			break;
		}
	}

	//Do the full error calculation if requested
	if ( ErrorMode > 0 )
	{
		//Do the full error calculation, either just for the variances (ErrorMode == 1) or for all covariances (ErrorMode == 2)
		bool justVariance = ( ErrorMode == 1 );
		fullErrors = new CovarianceMatrix( lastUnfoldingMatrix, inputSmearing, dataDistribution, unfoldedDistribution->Integral(), justVariance );

		//Read out the variance
		dagostiniVariance = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
		for ( unsigned int binIndex = 0; binIndex < indexCalculator->GetBinNumber(); binIndex++ )
		{
			dagostiniVariance[ binIndex ] = fullErrors->GetElement( binIndex, binIndex );
		}

		if ( justVariance )
		{
			delete fullErrors;
		}
	}
	delete lastUnfoldingMatrix;

	//Close the debug file
	if ( debug )
	{
		debugFile->Close();
	}
}

//Perform a closure test
//Unfold the MC reco distribution with the corresponding truth information as a prior
//It should give the truth information back exactly...
bool IterativeUnfolding::ClosureTest( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	//Use the truth distribution as the prior
	Distribution * priorDistribution = truthDistribution;

	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Make a pointer for the iteration result
	Distribution * unfoldedReconstructedDistribution;

	//Iterate
	UnfoldingMatrix * lastUnfoldingMatrix;
	for ( unsigned int iteration = 0; iteration < MostIterations; iteration++ )
	{
		//Smooth the prior distribution, is asked. Don't smooth the truth
		if ( WithSmoothing && iteration != 0 )
		{
			priorDistribution->Smooth();
		}

		//Unfold
		lastUnfoldingMatrix = new UnfoldingMatrix( inputSmearing, priorDistribution );
		unfoldedReconstructedDistribution = new Distribution( reconstructedDistribution, lastUnfoldingMatrix );
		delete lastUnfoldingMatrix;

		//Compare with previous distribution
		double chi2, kolmogorov;
		distributionComparison->CompareDistributions( unfoldedReconstructedDistribution, priorDistribution, chi2, kolmogorov, false );

		//Reset for next iteration
		if ( iteration != 0 )
		{
			//Don't delete the MC truth
			delete priorDistribution;
		}
		priorDistribution = unfoldedReconstructedDistribution;

		//Check for termination conditions
		if ( chi2 < ChiSquaredThreshold || kolmogorov > KolmogorovThreshold || iteration == MostIterations - 1 )
		{
			break;
		}
	}

	//Compare with truth distribution
	double chi2Reference, kolmogorovReference;
	distributionComparison->CompareDistributions( truthDistribution, unfoldedReconstructedDistribution, chi2Reference, kolmogorovReference, false, true );

	//Output result
	double binNumber = (double)indexCalculator->GetBinNumber();
	delete unfoldedReconstructedDistribution;
	if ( chi2Reference == 0.0 && kolmogorovReference == 1.0 )
	{
		cout << "Perfect closure test: chi squared = " << chi2Reference << " and K-S probability = " << kolmogorovReference << ". Nice one!" << endl;
		return true;
	}
	else if ( chi2Reference <  binNumber && kolmogorovReference > 1.0 / binNumber )
	{
		cout << "Closure test passed: chi squared = " << chi2Reference << " and K-S probability = " << kolmogorovReference << endl;
		return true;
	}
	else
	{
		cout << "Closure test failed: chi squared = " << chi2Reference << " and K-S probability = " << kolmogorovReference << endl;
		return false;
	}
}

//Perform an unfolding cross-check
//Use MC truth A as a prior to unfold MC reco B
//Iterations cease when result is sufficiently close to MC truth B (passed as argument)
//Returns the number of iterations required. Convergence criteria as output arguments
unsigned int IterativeUnfolding::MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	//Extrapolate the number of missed events in the data
	dataDistribution->SetBadBin( totalMissed / ( totalPaired + totalFake ) );

	//Use the truth distribution as the prior
	Distribution * priorDistribution = truthDistribution;

	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Compare the uncorrected data to the truth
	double lastChiSquared, lastKolmogorov;
	distributionComparison->CompareDistributions( truthDistribution, dataDistribution, lastChiSquared, lastKolmogorov, true );
	double chiSquaredResult = lastChiSquared;
	double kolmogorovResult = lastKolmogorov;
	cout << "------------- Cross-Check -------------" << endl;
	cout << "0: " << lastChiSquared << ", " << lastKolmogorov;

	//Iterate, making new distribution from data, old distribution and smearing matrix
	Distribution * adjustedDistribution;
	UnfoldingMatrix * lastUnfoldingMatrix;
	for ( unsigned int iteration = 0; iteration < MAX_ITERATIONS_FOR_CROSS_CHECK; iteration++ )
	{
		//Smooth the prior distribution, is asked. Don't smooth the truth
		if ( WithSmoothing && iteration != 0 )
		{
			priorDistribution->Smooth();
		}

		//Iterate
		lastUnfoldingMatrix = new UnfoldingMatrix( inputSmearing, priorDistribution );
		adjustedDistribution = new Distribution( dataDistribution, lastUnfoldingMatrix );
		delete lastUnfoldingMatrix;

		//Compare with reference distribution
		double referenceChi2, referenceKolmogorov;
		distributionComparison->CompareDistributions( adjustedDistribution, ReferenceDistribution, referenceChi2, referenceKolmogorov, true );


		//Compare with last iteration
		double chi2, kolmogorov;
		distributionComparison->CompareDistributions( adjustedDistribution, priorDistribution, chi2, kolmogorov, false );

		//Reset for next iteration
		if ( iteration != 0 )
		{
			//Don't delete the MC truth
			delete priorDistribution;
		}
		priorDistribution = adjustedDistribution;

		//Check to see if things have got worse
		if ( referenceChi2 > lastChiSquared || referenceKolmogorov < lastKolmogorov || iteration == MAX_ITERATIONS_FOR_CROSS_CHECK - 1 || ( referenceChi2 == lastChiSquared && referenceKolmogorov == lastKolmogorov ) )
		{
			//Return the criteria
			ChiSquaredThreshold = chiSquaredResult;
			KolmogorovThreshold = kolmogorovResult;
			cout << " <--" << endl << iteration + 1 << ": " << referenceChi2 << ", " << referenceKolmogorov << endl;
			cout << "-------------------------------------" << endl;
			delete adjustedDistribution;
			return iteration;
		}
		else
		{
			//Update the last values
			lastChiSquared = referenceChi2;
			lastKolmogorov = referenceKolmogorov;
			chiSquaredResult = chi2;
			kolmogorovResult = kolmogorov;
			cout << endl << iteration + 1 << ": " << referenceChi2 << ", " << referenceKolmogorov;
		}
	}

	//It should be impossible to get this far, but to prevent -Wall from complaining
	return MAX_ITERATIONS_FOR_CROSS_CHECK;
}

//Retrieve a TH1F* containing the unfolded data
//distribution, with or without errors
//NB: the error calculation is only performed
//when you run the method with errors for the first time
TH1F * IterativeUnfolding::GetUnfoldedHistogram( string Name, string Title, bool Normalise )
{
	return unfoldedDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Retrieve the smearing matrix used
TH2F * IterativeUnfolding::GetSmearingMatrix( string Name, string Title )
{
	return inputSmearing->MakeRootHistogram( Name, Title );
}

//Retrieve the truth distribution
TH1F * IterativeUnfolding::GetTruthHistogram( string Name, string Title, bool Normalise )
{
	return truthDistribution->MakeRootHistogram( Name, Title, Normalise );
}
Distribution * IterativeUnfolding::GetTruthDistribution()
{
	return truthDistribution;
}

//Retrieve the uncorrected data distribution
TH1F * IterativeUnfolding::GetUncorrectedDataHistogram( string Name, string Title, bool Normalise )
{
	return dataDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Handy for error calculation
vector< double > IterativeUnfolding::SumOfDataWeightSquares()
{
	return sumOfDataWeightSquares;
}
vector< double > IterativeUnfolding::DAgostiniVariance()
{
        return dagostiniVariance;
}
TH2F * IterativeUnfolding::DAgostiniCovariance( string Name, string Title )
{
	return fullErrors->MakeRootHistogram( Name, Title );
}
