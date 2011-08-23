/**
  @class BinByBinUnfolding

  Unfold the distribution with the simple bin-by-bin method

  @author Benjamin M Wynne bwynne@cern.ch
  @date 07-07-2011
 */

#include "BinByBinUnfolding.h"
#include <iostream>
#include <cstdlib>
#include <sstream>

//Default constructor - useless
BinByBinUnfolding::BinByBinUnfolding()
{
}

//N-Dimensional version
BinByBinUnfolding::BinByBinUnfolding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID )
{
	name = Name;
	uniqueID = UniqueID;
	totalPaired = 0.0;
	totalFake = 0.0;
	totalMissed = 0.0;
	isClone = false;

	indexCalculator = DistributionIndices;
	dataDistribution = new Distribution( indexCalculator );
	unfoldedDistribution = 0;
	truthDistribution = new Distribution( indexCalculator );
	reconstructedDistribution = new Distribution( indexCalculator );
	sumOfDataWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
	distributionComparison = new Comparison( Name, UniqueID );

	//Make the vectors for storing the bin-by-bin correction ratios (include bad bins)
	truthBinSums = vector< double >( indexCalculator->GetBinNumber() + 1, 0.0 );
	recoBinSums = vector< double >( indexCalculator->GetBinNumber() + 1, 0.0 );
}

//To be used with Clone
BinByBinUnfolding::BinByBinUnfolding( IIndexCalculator * DistributionIndices, string Name, unsigned int UniqueID,
		Comparison * SharedComparison, Distribution * SharedTruth, vector< double > & SharedTruthSums, vector< double > & SharedRecoSums, double PairedMC, double MissedMC, double FakeMC )
{
	name = Name;
	uniqueID = UniqueID;
	totalPaired = PairedMC;
	totalFake = FakeMC;
	totalMissed = MissedMC;
	isClone = true;

	indexCalculator = DistributionIndices;
	dataDistribution = new Distribution( indexCalculator );
	unfoldedDistribution = 0;
	truthDistribution = SharedTruth;
	reconstructedDistribution = new Distribution( indexCalculator );
	sumOfDataWeightSquares = vector< double >( indexCalculator->GetBinNumber(), 0.0 );
	distributionComparison = SharedComparison;

	//Make the vectors for storing the bin-by-bin correction ratios (include bad bins)
	truthBinSums = SharedTruthSums;
	recoBinSums = SharedRecoSums;;
}

//Destructor
BinByBinUnfolding::~BinByBinUnfolding()
{
	delete dataDistribution;
	delete reconstructedDistribution;
	if ( unfoldedDistribution )
	{
		delete unfoldedDistribution;
	}

	if ( isClone )
	{
		delete indexCalculator;
	}
	else
	{
		delete truthDistribution;
		delete distributionComparison;
	}

	sumOfDataWeightSquares.clear();
	truthBinSums.clear();
	recoBinSums.clear();
}

//Make another instance of the ICorrection which shares the smearing matrix
BinByBinUnfolding * BinByBinUnfolding::CloneShareSmearingMatrix()
{
	return new BinByBinUnfolding( indexCalculator->Clone(), name, uniqueID + 1, distributionComparison, truthDistribution, truthBinSums, recoBinSums, totalPaired, totalMissed, totalFake );
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void BinByBinUnfolding::StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreEvent( Truth, TruthWeight );
		reconstructedDistribution->StoreEvent( Reco, RecoWeight );
	}

	totalPaired += TruthWeight;
	truthBinSums[ indexCalculator->GetIndex( Truth ) ] += TruthWeight;
	recoBinSums[ indexCalculator->GetIndex( Reco ) ] += RecoWeight;
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void BinByBinUnfolding::StoreUnreconstructedTruth( vector< double > Truth, double Weight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreEvent( Truth, Weight );
		reconstructedDistribution->StoreBadEvent( Weight );
	}

	totalMissed += Weight;
	truthBinSums[ indexCalculator->GetIndex( Truth ) ] += Weight;
	recoBinSums[ recoBinSums.size() - 1 ] += Weight;
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void BinByBinUnfolding::StoreReconstructedFake( vector< double > Reco, double Weight, bool UseInPrior )
{
	if ( UseInPrior )
	{
		truthDistribution->StoreBadEvent( Weight );
		reconstructedDistribution->StoreEvent( Reco, Weight );
	}

	totalFake += Weight;
	truthBinSums[ truthBinSums.size() - 1 ] += Weight;
	recoBinSums[ indexCalculator->GetIndex( Reco ) ] += Weight;
}

//Store a value from the uncorrected data distribution
void BinByBinUnfolding::StoreDataValue( vector< double > Data, double Weight )
{
	dataDistribution->StoreEvent( Data, Weight );
	sumOfDataWeightSquares[ indexCalculator->GetIndex( Data ) ] += ( Weight * Weight );
}

//Once all data is stored, run the unfolding
//The arguments are all dummies
void BinByBinUnfolding::Correct( unsigned int MostIterations, unsigned int ErrorMode, bool WithSmoothing )
{
	//Extrapolate the number of missed events in the data
	dataDistribution->SetBadBin( totalMissed / ( totalPaired + totalFake ) );	

	//Work out the truth/reco ratios
	vector< double > binWeights( truthBinSums.size(), 0.0 );
	for ( unsigned int binIndex = 0; binIndex < truthBinSums.size(); binIndex++ )
	{
		if ( truthBinSums[ binIndex ] == 0.0 )
		{
			binWeights[ binIndex ] = 0.0;
		}
		else if ( recoBinSums[ binIndex ] == 0.0 )
		{
			binWeights[ binIndex ] = 1.0;
		}
		else
		{
			binWeights[ binIndex ] = truthBinSums[ binIndex ] / recoBinSums[ binIndex ];
		}
	}

	//Bin-by-bin
	unfoldedDistribution = new Distribution( dataDistribution, binWeights );
}

//Perform a closure test
//Unfold the MC reco distribution with the corresponding truth information as a prior
//It should give the truth information back exactly...
bool BinByBinUnfolding::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
{
	//Work out the truth/reco ratios
	vector< double > binWeights( truthBinSums.size(), 0.0 );
	for ( unsigned int binIndex = 0; binIndex < truthBinSums.size(); binIndex++ )
	{
		if ( truthBinSums[ binIndex ] == 0.0 )
		{
			binWeights[ binIndex ] = 0.0;
		}       
		else if ( recoBinSums[ binIndex ] == 0.0 )
		{
			binWeights[ binIndex ] = 1.0;
		}       
		else
		{
			binWeights[ binIndex ] = truthBinSums[ binIndex ] / recoBinSums[ binIndex ];
		}
	}

	//Bin-by-bin
	Distribution * unfoldedReconstructedDistribution = new Distribution( reconstructedDistribution, binWeights );

	//Compare with truth distribution
	double chi2Reference, kolmogorovReference;
	distributionComparison->CompareDistributions( truthDistribution, unfoldedReconstructedDistribution, chi2Reference, kolmogorovReference, true );

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
//Just a dummy since there is no iteration
unsigned int BinByBinUnfolding::MonteCarloCrossCheck( Distribution * ReferenceDistribution, bool WithSmoothing )
{
	return 1;
}

//Retrieve a TH1F* containing the correctted data distribution
TH1F * BinByBinUnfolding::GetCorrectedHistogram( string Name, string Title, bool Normalise )
{
	return unfoldedDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Retrieve the smearing matrix used
TH2F * BinByBinUnfolding::GetSmearingMatrix( string Name, string Title )
{
	//Make a matrix to store the bin-by-bin ratios
	TH2F * binByBinMatrix = new TH2F( Name.c_str(), Title.c_str(), truthBinSums.size(), 0.0, (double)truthBinSums.size(), truthBinSums.size(), 0.0, (double)truthBinSums.size() );

	//Work out the truth/reco ratios and store
	for ( unsigned int binIndex = 0; binIndex < truthBinSums.size(); binIndex++ )
	{
		int twoDimIndex = binByBinMatrix->GetBin( binIndex, binIndex, 0 );

		binByBinMatrix->SetBinContent( twoDimIndex, truthBinSums[ binIndex ] / recoBinSums[ binIndex ] );
	}

	return binByBinMatrix;
}

//Retrieve the truth distribution
TH1F * BinByBinUnfolding::GetTruthHistogram( string Name, string Title, bool Normalise )
{
	return truthDistribution->MakeRootHistogram( Name, Title, Normalise );
}
Distribution * BinByBinUnfolding::GetTruthDistribution()
{
	return truthDistribution;
}

//Retrieve the uncorrected data distribution
TH1F * BinByBinUnfolding::GetUncorrectedHistogram( string Name, string Title, bool Normalise )
{
	return dataDistribution->MakeRootHistogram( Name, Title, Normalise );
}

//Handy for error calculation
vector< double > BinByBinUnfolding::Variances()
{
	return sumOfDataWeightSquares;
}
TH2F * BinByBinUnfolding::DAgostiniCovariance( string Name, string Title )
{
	return 0;
}
