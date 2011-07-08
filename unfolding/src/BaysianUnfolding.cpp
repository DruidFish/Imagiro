/**
  @class BayesianUnfolding

  The class that manages unfolding a distribution, and tests the quality of the unfolding

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "BayesianUnfolding.h"
#include "UniformIndices.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <sstream>

const unsigned int MAX_ITERATIONS_FOR_CROSS_CHECK = 10;

//Default constructor - useless
BayesianUnfolding::BayesianUnfolding()
{
}

//N-Dimensional version
BayesianUnfolding::BayesianUnfolding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID )
{
	name = Name;
	uniqueID = UniqueID;
	totalPaired = 0.0;
	totalFake = 0.0;
	totalMissed = 0.0;

	indexCalculator = DistributionIndices;
	inputSmearing = new SmearingMatrix( indexCalculator );
	dataDistribution = new Distribution( indexCalculator );
	unfoldedDistribution = new Distribution( indexCalculator );
	truthDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = new Distribution( indexCalculator );
	sumOfDataWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
	distributionComparison = new Comparison( Name, UniqueID );
}

//Destructor
BayesianUnfolding::~BayesianUnfolding()
{
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
void BayesianUnfolding::StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
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
void BayesianUnfolding::StoreUnreconstructedTruth( vector< double > Truth, double Weight, bool UseInPrior )
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
void BayesianUnfolding::StoreReconstructedFake( vector< double > Reco, double Weight, bool UseInPrior )
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
void BayesianUnfolding::StoreDataValue( vector< double > Data, double Weight )
{
	dataDistribution->StoreEvent( Data, Weight );
	sumOfDataWeightSquares[ indexCalculator->GetIndex( Data ) ] += ( Weight * Weight );
}

//Once all data is stored, run the unfolding
//You can specify when the iterations should end,
//with an upper limit on iteration number
void BayesianUnfolding::Correct( unsigned int MostIterations, unsigned int ErrorMode, bool WithSmoothing )
{
	//Extrapolate the number of missed events in the data
	dataDistribution->SetBadBin( totalMissed / ( totalPaired + totalFake ) );

	//Use the truth distribution as the prior
	Distribution * priorDistribution = truthDistribution;

	//Finalise the smearing matrix
	inputSmearing->Finalise();

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

		//Reset for next iteration
		if ( iteration != 0 )
		{
			//Don't delete the MC truth
			delete priorDistribution;
		}
		priorDistribution = unfoldedDistribution;
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
}

//Perform a closure test
//Unfold the MC reco distribution with the corresponding truth information as a prior
//It should give the truth information back exactly...
bool BayesianUnfolding::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
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

		//Reset for next iteration
		if ( iteration != 0 )
		{
			//Don't delete the MC truth
			delete priorDistribution;
		}
		priorDistribution = unfoldedReconstructedDistribution;
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
unsigned int BayesianUnfolding::MonteCarloCrossCheck( Distribution * ReferenceDistribution, bool WithSmoothing )
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
			cout << endl << iteration + 1 << ": " << referenceChi2 << ", " << referenceKolmogorov;
		}
	}

	//It should be impossible to get this far, but to prevent -Wall from complaining
	return MAX_ITERATIONS_FOR_CROSS_CHECK;
}

//Retrieve a TH1F* containing the correctted data distribution
TH1F * BayesianUnfolding::GetCorrectedHistogram( string Name, string Title, bool Normalise )
{
	return unfoldedDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Retrieve the smearing matrix used
TH2F * BayesianUnfolding::GetSmearingMatrix( string Name, string Title )
{
	return inputSmearing->MakeRootHistogram( Name, Title );
}

//Retrieve the truth distribution
TH1F * BayesianUnfolding::GetTruthHistogram( string Name, string Title, bool Normalise )
{
	return truthDistribution->MakeRootHistogram( Name, Title, Normalise );
}
Distribution * BayesianUnfolding::GetTruthDistribution()
{
	return truthDistribution;
}

//Retrieve the uncorrected data distribution
TH1F * BayesianUnfolding::GetUncorrectedHistogram( string Name, string Title, bool Normalise )
{
	return dataDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Handy for error calculation
vector< double > BayesianUnfolding::Variances()
{
	if ( dagostiniVariance.size() > 0 )
	{
        	return dagostiniVariance;
	}
	else
	{
		return sumOfDataWeightSquares;
	}
}
TH2F * BayesianUnfolding::DAgostiniCovariance( string Name, string Title )
{
	return fullErrors->MakeRootHistogram( Name, Title );
}
