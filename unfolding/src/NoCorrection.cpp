/**
  @class NoCorrection

  Just creates distributions from the given data, with no attempt at correction

  @author Benjamin M Wynne bwynne@cern.ch
  @date 07-07-2011
 */

#include "NoCorrection.h"
#include <iostream>
#include <cstdlib>

//Default constructor - useless
NoCorrection::NoCorrection()
{
}

//Constructor taking an IIndexCalculator to define the bins
NoCorrection::NoCorrection( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID )
{
	//Initialisations
	name = Name;
	uniqueID = UniqueID;
	totalPaired = 0.0;
	totalFake = 0.0;
	totalMissed = 0.0;
	isClone = false;

	//Make new vectors just with the one entry
	indexCalculator = DistributionIndices;
	inputSmearing = new SmearingMatrix( indexCalculator );
	inputDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = new Distribution( indexCalculator );
	sumOfInputWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
}

//For use with Clone
NoCorrection::NoCorrection( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID,
	       Distribution * SharedReconstructed, SmearingMatrix * SharedSmearing, double PairedMC, double FakeMC, double MissedMC )
{
	//Initialisations
	name = Name;
	uniqueID = UniqueID;
	totalPaired = PairedMC;
	totalFake = FakeMC;
	totalMissed = MissedMC;
	isClone = true;

	//Make new vectors just with the one entry
	indexCalculator = DistributionIndices;
	inputSmearing = SharedSmearing;
	inputDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = SharedReconstructed;
	sumOfInputWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
}


//Destructor
NoCorrection::~NoCorrection()
{
	delete inputDistribution;
	if ( isClone )
	{
		delete indexCalculator;
	}
	else
	{
		delete reconstructedDistribution;
		delete inputSmearing;
	}
}

//Make another instance of the ICorrection which sihares the smearing matrix
NoCorrection * NoCorrection::CloneShareSmearingMatrix()
{
	return new NoCorrection( indexCalculator, name, uniqueID + 1, reconstructedDistribution, inputSmearing, totalPaired, totalFake, totalMissed );
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void NoCorrection::StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if (UseInPrior)
	{
		reconstructedDistribution->StoreEvent( Reco, RecoWeight );
	}

	totalPaired += TruthWeight;
	inputSmearing->StoreTruthRecoPair( Truth, Reco, TruthWeight, RecoWeight );
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void NoCorrection::StoreUnreconstructedTruth( vector< double > Truth, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		reconstructedDistribution->StoreBadEvent( Weight );
	}

	totalMissed += Weight;
	inputSmearing->StoreUnreconstructedTruth( Truth, Weight );
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void NoCorrection::StoreReconstructedFake( vector< double > Reco, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		reconstructedDistribution->StoreEvent( Reco, Weight );
	}

	totalFake += Weight;
	inputSmearing->StoreReconstructedFake( Reco, Weight );
}

//Store a value from the uncorrected data distribution
void NoCorrection::StoreDataValue( vector< double > Input, double Weight )
{
	inputDistribution->StoreEvent( Input, Weight );
	sumOfInputWeightSquares[ indexCalculator->GetIndex( Input ) ] += ( Weight * Weight );
}

//Dummy, since nothing is happening
void NoCorrection::Correct( unsigned int MostIterations, unsigned int ErrorMode, bool WithSmoothing )
{
	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Extrapolate the number of missed events in the data
	inputDistribution->SetBadBin( totalMissed / ( totalPaired + totalFake ) );	
}

//Perform a closure test
//Dummy, since nothing is happening
bool NoCorrection::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
{
	cout << "No correction, therefore no closure test to perform" << endl;
	return true;
}

//Perform an unfolding cross-check
//Dummy, since folding is not iterative
unsigned int NoCorrection::MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing )
{
	return 1;
}

//Retrieve the uncorrected distribution
TH1F * NoCorrection::GetCorrectedHistogram( string Name, string Title, bool Normalise )
{
	return inputDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Retrieve the smearing matrix used
TH2F * NoCorrection::GetSmearingHistogram( string Name, string Title )
{
	return inputSmearing->MakeRootHistogram( Name, Title );
}
SmearingMatrix * NoCorrection::GetSmearingMatrix()
{
	return inputSmearing;
}

//Retrieve the reconstructed distribution
TH1F * NoCorrection::GetTruthHistogram( string Name, string Title, bool Normalise )
{
	return reconstructedDistribution->MakeRootHistogram( Name, Title, Normalise );
}
Distribution * NoCorrection::GetTruthDistribution()
{
	return reconstructedDistribution;
}

//Retrieve the uncorrected data distribution
TH1F * NoCorrection::GetUncorrectedHistogram( string Name, string Title, bool Normalise )
{
	return inputDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Handy for error calculation
vector< double > NoCorrection::Variances()
{
	return sumOfInputWeightSquares;
}
TH2F * NoCorrection::DAgostiniCovariance( string Name, string Title )
{
	return 0;
}
