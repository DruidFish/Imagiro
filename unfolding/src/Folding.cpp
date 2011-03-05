/**
  @class Folding

  The class that applies a smearing matrix to a "truth" distribution, to simulate the detector effects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 11-02-2011
 */

#include "Folding.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>

//Default constructor - useless
Folding::Folding()
{
}

//Constructor taking the required bin number,
//minimum and maximum of the output distribution as arguments
//NB: the unfolding scales roughly with bin number ^ 2, the
//error calculation scales roughly with bin number ^ 3.
Folding::Folding( int BinNumber, double Minimum, double Maximum, string Name, int UniqueID ) : name(Name), uniqueID(UniqueID)
{
	//Make new vectors just with the one entry
	indexCalculator = new Indices( vector<int>( 1, BinNumber ), vector<double>( 1, Minimum ), vector<double>( 1, Maximum ) );
	inputSmearing = new SmearingMatrix(indexCalculator);
	truthDistribution = new Distribution(indexCalculator);
	inputDistribution = new Distribution(indexCalculator);
	reconstructedDistribution = new Distribution(indexCalculator);
	sumOfInputWeightSquares = vector<double>( indexCalculator->GetBinNumber(), 0.0 );
	distributionComparison = new Comparison( Name, UniqueID );
}

//N-Dimensional version
Folding::Folding( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima, string Name, int UniqueID ) : name(Name), uniqueID(UniqueID)
{
	indexCalculator = new Indices( BinNumbers, Minima, Maxima );
	inputSmearing = new SmearingMatrix(indexCalculator);
	truthDistribution = new Distribution(indexCalculator);
	inputDistribution = new Distribution(indexCalculator);
	reconstructedDistribution = new Distribution(indexCalculator);
	sumOfInputWeightSquares = vector<double>( indexCalculator->GetBinNumber(), 0.0 );
	distributionComparison = new Comparison( Name, UniqueID );
}

//Destructor
Folding::~Folding()
{
	delete indexCalculator;
	delete inputSmearing;
	delete truthDistribution;
	delete inputDistribution;
	delete reconstructedDistribution;
	delete smearedDistribution;
	delete distributionComparison;
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void Folding::StoreTruthRecoPair( double Truth, double Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreEvent( vector<double>( 1, Truth ), TruthWeight );
		reconstructedDistribution->StoreEvent( vector<double>( 1, Reco ), RecoWeight );
	}
	inputSmearing->StoreTruthRecoPair( vector<double>( 1, Truth ), vector<double>( 1, Reco ), TruthWeight, RecoWeight );
}

//N-Dimensional version
void Folding::StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double TruthWeight, double RecoWeight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreEvent( Truth, TruthWeight );
		reconstructedDistribution->StoreEvent( Reco, RecoWeight );
	}
	inputSmearing->StoreTruthRecoPair( Truth, Reco, TruthWeight, RecoWeight );
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void Folding::StoreUnreconstructedTruth( double Truth, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreEvent( vector<double>( 1, Truth ), Weight );
	}
	inputSmearing->StoreUnreconstructedTruth( vector<double>( 1, Truth ), Weight );
}

//N-Dimensional version
void Folding::StoreUnreconstructedTruth( vector<double> Truth, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		truthDistribution->StoreEvent( Truth, Weight );
	}
	inputSmearing->StoreUnreconstructedTruth( Truth, Weight );
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void Folding::StoreReconstructedFake( double Reco, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		reconstructedDistribution->StoreEvent( vector<double>( 1, Reco ), Weight );
	}
	inputSmearing->StoreReconstructedFake( vector<double>( 1, Reco ), Weight );
}

//N-Dimensional version
void Folding::StoreReconstructedFake( vector<double> Reco, double Weight, bool UseInPrior )
{
	if (UseInPrior)
	{
		reconstructedDistribution->StoreEvent( Reco, Weight );
	}
	inputSmearing->StoreReconstructedFake( Reco, Weight );
}

//Store a value from the uncorrected data distribution
void Folding::StoreValueToFold( double Input, double Weight )
{
	vector<double> inputVector( 1, Input );
	inputDistribution->StoreEvent( inputVector, Weight );
	sumOfInputWeightSquares[ indexCalculator->GetIndex( inputVector ) ] += Weight;
}

//N-Dimensional version
void Folding::StoreValueToFold( vector<double> Input, double Weight )
{
	inputDistribution->StoreEvent( Input, Weight );
	sumOfInputWeightSquares[ indexCalculator->GetIndex( Input ) ] += Weight;
}

//Smear the input distribution
void Folding::Fold()
{
	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Make the smeared distribution
	smearedDistribution = new Distribution( inputDistribution, inputSmearing );
}

//Perform a closure test
//Fold the MC truth information - should give the MC reco exactly
void Folding::ClosureTest()
{
	//Finalise the smearing matrix
	inputSmearing->Finalise();

	//Fold the truth distribution
	Distribution * smearedTruthDistribution = new Distribution( truthDistribution, inputSmearing );

	//Compare with reconstructed distribution
	double chi2, kolmogorov;
	distributionComparison->CompareDistributions( reconstructedDistribution, smearedTruthDistribution, chi2, kolmogorov, false );

	//Output result
	if ( chi2 < 1.0 && kolmogorov > 0.9 )
	{
		cout << "Closure test passed: chi squared = " << chi2 << " and K-S probability = " << kolmogorov << endl;
	}       
	else
	{
		cout << "Closure test failed: chi squared = " << chi2 << " and K-S probability = " << kolmogorov << endl;
	}
}

//Retrieve the folded distribution
TH1F * Folding::GetFoldedHistogram( string Name, string Title, bool WithErrors )
{
	return smearedDistribution->MakeRootHistogram( Name, Title );
}

//Retrieve the smearing matrix used
TH2F * Folding::GetSmearingMatrix( string Name, string Title )
{
	return inputSmearing->MakeRootHistogram( Name, Title );
}

//Retrieve the truth distribution
TH1F * Folding::GetReconstructedHistogram( string Name, string Title )
{
	return reconstructedDistribution->MakeRootHistogram( Name, Title );
}

//Retrieve the truth distribution
TH1F * Folding::GetTruthHistogram( string Name, string Title )
{
	return truthDistribution->MakeRootHistogram( Name, Title );
}
Distribution * Folding::GetTruthDistribution()
{
	return truthDistribution;
}

//Retrieve the uncorrected data distribution
TH1F * Folding::GetInputHistogram( string Name, string Title )
{
	return inputDistribution->MakeRootHistogram( Name, Title );
}

//Handy for error calculation
vector<double> Folding::SumOfInputWeightSquares()
{
	return sumOfInputWeightSquares;
}
