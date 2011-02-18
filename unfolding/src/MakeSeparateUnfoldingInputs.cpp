/**
  @class MakeSeparateUnfoldingInputs

  OUTDATED
  The class to make the inputs for an unfolding process in situations where only some fraction of the events are available

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#include "MakeSeparateUnfoldingInputs.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>

//Default constructor - useless
MakeSeparateUnfoldingInputs::MakeSeparateUnfoldingInputs()
{
}

//Constructor taking the required bin number,
//minimum and maximum of the output distribution as arguments
//NB: the unfolding scales roughly with bin number ^ 2, the
//error calculation scales roughly with bin number ^ 3.
MakeSeparateUnfoldingInputs::MakeSeparateUnfoldingInputs( int BinNumber, double Minimum, double Maximum )
{
	indexCalculator = new Indices( vector<int>( 1, BinNumber ), vector<double>( 1, Minimum ), vector<double>( 1, Maximum ) );
	inputSmearing = new SmearingMatrix(indexCalculator);
	truth = new Distribution(indexCalculator);
	data = new Distribution(indexCalculator);
}

//N-Dimensional version
MakeSeparateUnfoldingInputs::MakeSeparateUnfoldingInputs( vector<int> BinNumbers, vector<double> Minima, vector<double> Maxima )
{
	indexCalculator = new Indices( BinNumbers, Minima, Maxima );
	inputSmearing = new SmearingMatrix(indexCalculator);
        truth = new Distribution(indexCalculator);
        data = new Distribution(indexCalculator);
}

//Destructor
MakeSeparateUnfoldingInputs::~MakeSeparateUnfoldingInputs()
{
	delete indexCalculator;
	delete inputSmearing;
	delete truth;
	delete data;
}

//Use this method to supply a value from the truth
//distribution, and the corresponding reconstructed
//value
//NB: These values must both come from the SAME
//Monte Carlo event, or the whole process is meaningless
void MakeSeparateUnfoldingInputs::StoreTruthRecoPair( double Truth, double Reco, double Weight )
{
	truth->StoreEvent( vector<double>( 1, Truth ), Weight );
	inputSmearing->StoreTruthRecoPair( vector<double>( 1, Truth ), vector<double>( 1, Reco ), Weight );
}

//N-Dimensional version
void MakeSeparateUnfoldingInputs::StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double Weight )
{
	truth->StoreEvent( Truth, Weight );
        inputSmearing->StoreTruthRecoPair( Truth, Reco, Weight );
}

//If an MC event is not reconstructed at all, use this
//method to store the truth value alone
void MakeSeparateUnfoldingInputs::StoreUnreconstructedTruth( double Truth, double Weight )
{
	truth->StoreEvent( vector<double>( 1, Truth ), Weight );
	inputSmearing->StoreUnreconstructedTruth( vector<double>( 1, Truth ), Weight );
}

//N-Dimensional version
void MakeSeparateUnfoldingInputs::StoreUnreconstructedTruth( vector<double> Truth, double Weight )
{
	truth->StoreEvent( Truth, Weight );
        inputSmearing->StoreUnreconstructedTruth( Truth, Weight );
}

//If there is a fake reconstructed event with no
//corresponding truth, use this method
void MakeSeparateUnfoldingInputs::StoreReconstructedFake( double Reco, double Weight )
{
        inputSmearing->StoreReconstructedFake( vector<double>( 1, Reco ), Weight );
}

//N-Dimensional version
void MakeSeparateUnfoldingInputs::StoreReconstructedFake( vector<double> Reco, double Weight )
{
	inputSmearing->StoreReconstructedFake( Reco, Weight );
}

//Store a value from the uncorrected data distribution
void MakeSeparateUnfoldingInputs::StoreDataValue( double Data, double Weight )
{
	data->StoreEvent( vector<double>( 1, Data ), Weight );
}

//N-Dimensional version
void MakeSeparateUnfoldingInputs::StoreDataValue( vector<double> Data, double Weight )
{
	data->StoreEvent( Data, Weight );
}

//Retrieve an (un-normalised) smearing matrix to use
//in a later calculation
TH2F * MakeSeparateUnfoldingInputs::UnnormalisedSmearingMatrix()
{
	return inputSmearing->MakeRootHistogram( "smearing", "Unnormalised smearing matrix", false );
}

//Retrieve an (un-normalised) truth distribution
//to use in a later calculation
TH1F * MakeSeparateUnfoldingInputs::UnnormalisedTruthDistribution()
{
	return truth->MakeRootHistogram( "truth", "Unnormalised truth distribution", false );
}

//Retrieve an (un-normalised) data distribution
//to use in a later calculation
TH1F * MakeSeparateUnfoldingInputs::UnnormalisedDataDistribution()
{
	return data->MakeRootHistogram( "data", "Unnormalised data distribution", false );
}
