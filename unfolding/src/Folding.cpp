/**
  @class Folding

  The class that applies a smearing matrix to a "truth" distribution, to simulate the detector effects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-02-2011
 */

#include "Folding.h"
#include <iostream>
#include <cstdlib>

//Default constructor - useless
Folding::Folding()
{
}

//Constructor taking an IIndexCalculator to define the bins
Folding::Folding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID )
{
	//Initialisations
	name = Name;
	uniqueID = UniqueID;
	totalPaired = 0.0;
	totalMissed = 0.0;
	totalFake = 0.0;

	//Make new vectors just with the one entry
	indexCalculator = DistributionIndices;
	inputSmearing = new SmearingMatrix( indexCalculator );
	truthDistribution = new Distribution( indexCalculator );
	inputDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = new Distribution( indexCalculator );
	sumOfInputWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
	distributionComparison = new Comparison( Name, UniqueID );
	smearedDistribution = 0;
}

//Destructor
Folding::~Folding()
{
	delete inputSmearing;
	delete truthDistribution;
	delete inputDistribution;
	delete reconstructedDistribution;
	if ( smearedDistribution )
	{
		delete smearedDistribution;
	}
	delete distributionComparison;
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void Folding::StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreEvent( Truth, TruthWeight );
		reconstructedDistribution->StoreEvent( Reco, RecoWeight );
	}

	totalPaired += TruthWeight;
	inputSmearing->StoreTruthRecoPair( Truth, Reco, TruthWeight, RecoWeight );
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void Folding::StoreUnreconstructedTruth( vector< double > Truth, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreEvent( Truth, Weight );
		reconstructedDistribution->StoreBadEvent( Weight );
	}

	totalMissed += Weight;
	inputSmearing->StoreUnreconstructedTruth( Truth, Weight );
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void Folding::StoreReconstructedFake( vector< double > Reco, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreBadEvent( Weight );
		reconstructedDistribution->StoreEvent( Reco, Weight );
	}

	totalFake += Weight;
	inputSmearing->StoreReconstructedFake( Reco, Weight );
}

//Store a value from the uncorrected data distribution
void Folding::StoreDataValue( vector< double > Input, double Weight )
{
	inputDistribution->StoreEvent( Input, Weight );
	sumOfInputWeightSquares[ indexCalculator->GetIndex( Input ) ] += ( Weight * Weight );
}

//Smear the input distribution
void Folding::Correct( unsigned int MostIterations, unsigned int ErrorMode, bool WithSmoothing )
{
	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Extrapolate the number of fake events in the input distribution
	inputDistribution->SetBadBin( totalFake / ( totalPaired + totalMissed ) );

	//Make the smeared distribution
	smearedDistribution = new Distribution( inputDistribution, inputSmearing );
}

//Perform a closure test
//Fold the MC truth information - should give the MC reco exactly
bool Folding::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
{
	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Fold the truth distribution
	Distribution * smearedTruthDistribution = new Distribution( truthDistribution, inputSmearing );

	//Compare with reconstructed distribution
	double chi2, kolmogorov;
	distributionComparison->CompareDistributions( reconstructedDistribution, smearedTruthDistribution, chi2, kolmogorov, false, true );

	//Output result
	cout << "Number of bins: " << indexCalculator->GetBinNumber() << endl;
	if ( chi2 == 0.0 && kolmogorov == 1.0 )
	{
		cout << "Perfect closure test: chi squared = " << chi2 << " and K-S probability = " << kolmogorov << ". Nice one!" << endl;
		return true;
	}
	else if ( chi2 < 10.0 && kolmogorov > 0.1 )
	{
		cout << "Closure test passed: chi squared = " << chi2 << " and K-S probability = " << kolmogorov << endl;
		return true;
	}
	else
	{
		cout << "Closure test failed: chi squared = " << chi2 << " and K-S probability = " << kolmogorov << endl;
		return false;
	}
}

//Perform an unfolding cross-check
//Dummy, since folding is not iterative
unsigned int Folding::MonteCarloCrossCheck( Distribution * ReferenceDistribution, bool WithSmoothing )
{
	return 1;
}

//Retrieve the folded distribution
TH1F * Folding::GetCorrectedHistogram( string Name, string Title, bool Normalise )
{
	return smearedDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Retrieve the smearing matrix used
TH2F * Folding::GetSmearingMatrix( string Name, string Title )
{
	return inputSmearing->MakeRootHistogram( Name, Title );
}

//Retrieve the reconstructed distribution
TH1F * Folding::GetTruthHistogram( string Name, string Title, bool Normalise )
{
	return reconstructedDistribution->MakeRootHistogram( Name, Title, Normalise );
}
Distribution * Folding::GetTruthDistribution()
{
	return reconstructedDistribution;
}

//Retrieve the uncorrected data distribution
TH1F * Folding::GetUncorrectedHistogram( string Name, string Title, bool Normalise )
{
	return inputDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Handy for error calculation
vector< double > Folding::Variances()
{
	return sumOfInputWeightSquares;
}
TH2F * Folding::DAgostiniCovariance( string Name, string Title )
{
	return 0;
}
