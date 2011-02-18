#include "XvsYNormalisedFolding.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "TFile.h"

using namespace std;

//Default constructor - useless
XvsYNormalisedFolding::XvsYNormalisedFolding()
{
}

//Constructor with the names to use for the variables
XvsYNormalisedFolding::XvsYNormalisedFolding( string XVariableName, string YVariableName, string PriorName,
		int XBinNumber, double XMinimum, double XMaximum,
		int YBinNumber, double YMinimum, double YMaximum,
		double ScaleFactor, int UniqueID ) : xName( XVariableName ), yName( YVariableName ), priorName( PriorName ), finalised( false ), uniqueID(UniqueID), scaleFactor(ScaleFactor)
{
	vector<double> minima, maxima;
	vector<int> binNumbers;

	//Store the x range
	minima.push_back( XMinimum );
	maxima.push_back( XMaximum );
	binNumbers.push_back( XBinNumber );

	//Make the x unfolder
	XFolder = new Folding( binNumbers, minima, maxima, xName + priorName, uniqueID );

	//Store the y range
	minima.push_back( YMinimum );
	maxima.push_back( YMaximum );
	binNumbers.push_back( YBinNumber );

	//Make the x vs y unfolder
	XvsYFolder = new Folding( binNumbers, minima, maxima, xName + "vs" + yName + priorName, uniqueID );

	//Set up the indices for the distributions
	DistributionIndices = new DataIndices( binNumbers, minima, maxima );
}

//Destructor
XvsYNormalisedFolding::~XvsYNormalisedFolding()
{
	delete DistributionIndices;
	delete XFolder;
	delete XvsYFolder;
	delete foldedDistribution;
	delete inputDistribution;
	delete reconstructedDistribution;
	delete smearingMatrix;
}

//Copy the object
IPlotMaker * XvsYNormalisedFolding::Clone( string NewPriorName )
{
	return new XvsYNormalisedFolding( xName, yName, NewPriorName,
			DistributionIndices->GetBinNumber(0) - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0],
			DistributionIndices->GetBinNumber(1) - 2, DistributionIndices->GetMinima()[1], DistributionIndices->GetMaxima()[1],
			scaleFactor, uniqueID );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XvsYNormalisedFolding::StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
	else
	{
		vector<double> truthValues, reconstructedValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xTruthValue = TruthInput->GetValue( xName );
		double xReconstructedValue = ReconstructedInput->GetValue( xName );
		double yTruthValue = TruthInput->GetValue( yName );
		double yReconstructedValue = ReconstructedInput->GetValue( yName );
		double truthWeight = TruthInput->EventWeight();
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x values
		truthValues.push_back( xTruthValue );
		reconstructedValues.push_back( xReconstructedValue );
		XFolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		//Store the y values
		truthValues.push_back( yTruthValue );
		reconstructedValues.push_back( yReconstructedValue );
		XvsYFolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );
	}
}
void XvsYNormalisedFolding::StoreMiss( InputNtuple * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
	else
	{
		vector<double> truthValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xTruthValue = TruthInput->GetValue( xName );
		double yTruthValue = TruthInput->GetValue( yName );
		double truthWeight = TruthInput->EventWeight();

		//Store the x value
		truthValues.push_back( xTruthValue );
		XFolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );

		//Store the y value
		truthValues.push_back( yTruthValue );
		XvsYFolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );
	}
}
void XvsYNormalisedFolding::StoreFake( InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XvsYNormalisedFolding" << endl;
		exit(1);
	}       
	else
	{
		vector<double> reconstructedValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( ReconstructedInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xReconstructedValue = ReconstructedInput->GetValue( xName );
		double yReconstructedValue = ReconstructedInput->GetValue( yName );
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x value
		reconstructedValues.push_back( xReconstructedValue );
		XFolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );

		//Store the y value
		reconstructedValues.push_back( yReconstructedValue );
		XvsYFolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );
	}
}
void XvsYNormalisedFolding::StoreData( InputNtuple * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XvsYNormalisedFolding" << endl;
		exit(1);
	}       
	else
	{
		vector<double> dataValues;

		//Retrieve the values from the Ntuple
		double xDataValue = DataInput->GetValue( xName );
		double yDataValue = DataInput->GetValue( yName );
		double dataWeight = DataInput->EventWeight();

		//Store the x value
		dataValues.push_back( xDataValue );
		XFolder->StoreValueToFold( dataValues, dataWeight );

		//Store the y value
		dataValues.push_back( yDataValue );
		XvsYFolder->StoreValueToFold( dataValues, dataWeight );

		//Store values for performing the delinearisation
		DistributionIndices->StoreDataValue( dataValues, dataWeight );
	} 
}

//Do the unfolding
void XvsYNormalisedFolding::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XvsYNormalisedFolding is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Closure tests
		XFolder->ClosureTest();
		XvsYFolder->ClosureTest();

		//Unfold the distributions
		XFolder->Fold();
		XvsYFolder->Fold();

		//Make some plot titles
		char uniqueIDString[10];
		sprintf( uniqueIDString, "%d", uniqueID );
		string XFullName = xName + priorName + uniqueIDString;
		string XFullTitle = xName + " using " + priorName;
		string XvsYName = xName + "vs" + yName + priorName + uniqueIDString;
		string XvsYTitle = xName + " vs " + yName + " using " + priorName;

		//Retrieve the results
		TH1F * XSmeared = XFolder->GetFoldedHistogram( XFullName + "Smeared", XFullTitle + " Smeared Distribution" );
		TH1F * XvsYSmeared = XvsYFolder->GetFoldedHistogram( XvsYName + "Smeared", XvsYTitle + " Smeared Distribution" );

		//Retrieve some other bits for debug
		TH1F * XNotSmeared = XFolder->GetInputHistogram( XFullName + "Input", XFullTitle + " Input Distribution" );
		TH1F * XvsYNotSmeared = XvsYFolder->GetInputHistogram( XvsYName + "Input", XvsYTitle + " Input Distribution" );

		TH1F * XTruth = XFolder->GetReconstructedHistogram( XFullName + "Reconstructed", XFullTitle + " Reconstructed Distribution" );
		TH1F * XvsYTruth = XvsYFolder->GetReconstructedHistogram( XvsYName + "Reconstructed", XvsYTitle + " Reconstructed Distribution" );

		//TH2F * XSmearing = XFolder->GetSmearingMatrix( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );
		TH2F * XvsYSmearing = XvsYFolder->GetSmearingMatrix( XvsYName + "Smearing", XvsYTitle + " Smearing Matrix" );

		//De-linearise the x vs y distributions
		TH1F * DelinearisedXvsYSmeared = Delinearise(XvsYSmeared);
		TH1F * DelinearisedXvsYNotSmeared = Delinearise(XvsYNotSmeared);
		TH1F * DelinearisedXvsYTruth = Delinearise(XvsYTruth);

		//Get the error vectors
		vector<double> XErrors = XFolder->SumOfInputWeightSquares();
		vector<double> XvsYErrors = XvsYFolder->SumOfInputWeightSquares();

		//Delinearise the x vs y errors
		vector<double> delinearisedXvsYErrors = DelineariseErrors( XvsYErrors );

		//Combine errors
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			//Add the errors from the x vs y distribution and the divisor
			//The formula is stright from ROOT::TH1::Divide, so I hope it's right
			double XBinValue = XSmeared->GetBinContent(binIndex);
			double XvsYBinValue = DelinearisedXvsYSmeared->GetBinContent(binIndex);
			double componentOne = XErrors[binIndex] * XvsYBinValue * XvsYBinValue;
			double componentTwo = delinearisedXvsYErrors[binIndex] * XBinValue * XBinValue;
			double componentThree = XBinValue * XBinValue;
			double combinedError = ( sqrt( componentOne + componentTwo ) / componentThree );
			combinedError *= scaleFactor;
			correctedInputErrors.push_back( combinedError );
		}

		//Normalise the x vs y distributions and scale appropriately
		DelinearisedXvsYSmeared->Divide( DelinearisedXvsYSmeared, XSmeared, scaleFactor, 1.0 );
		DelinearisedXvsYNotSmeared->Divide( DelinearisedXvsYNotSmeared, XNotSmeared, scaleFactor, 1.0 );
		DelinearisedXvsYTruth->Divide( DelinearisedXvsYTruth, XTruth, scaleFactor, 1.0 );

		//Free some memory
		delete XSmeared;
		delete XNotSmeared;
		delete XTruth;

		//Get the y range to plot
		double yMinimum = DistributionIndices->GetMinima()[1] * scaleFactor;
		double yMaximum = DistributionIndices->GetMaxima()[1] * scaleFactor;

		//Format and save the corrected distribution
		foldedDistribution = DelinearisedXvsYSmeared;
		foldedDistribution->SetXTitle( xName.c_str() );
		foldedDistribution->SetYTitle( yName.c_str() );
		foldedDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the uncorrected distribution
		inputDistribution = DelinearisedXvsYNotSmeared;
		inputDistribution->SetXTitle( xName.c_str() );
		inputDistribution->SetYTitle( yName.c_str() );
		inputDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the truth distribution
		reconstructedDistribution = DelinearisedXvsYTruth;
		reconstructedDistribution->SetXTitle( xName.c_str() );
		reconstructedDistribution->SetYTitle( yName.c_str() );
		reconstructedDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the smearing matrix
		smearingMatrix = XvsYSmearing;
		string smearingXTitle = xName + " vs " + yName + " Truth";
		string smearingYTitle = xName + " vs " + yName + " Reconstructed";
		smearingMatrix->SetXTitle( smearingXTitle.c_str() );
		smearingMatrix->SetYTitle( smearingYTitle.c_str() );

		//Bin-by-bin scaling of errors using the corrected data
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double errorScaleFactor = foldedDistribution->GetBinContent(binIndex) / inputDistribution->GetBinContent(binIndex);
			correctedInputErrors[binIndex] *= errorScaleFactor;
		}

		//Mark as done
		finalised = true;
	}
}

//Make a cross-check with MC
int XvsYNormalisedFolding::MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	//Meaningless values
	ChiSquaredThreshold = 10.0;
	KolmogorovThreshold = 0.1;
	return 10;
}

//Return some plots
TH1F * XvsYNormalisedFolding::CorrectedHistogram()
{
	if ( finalised )
	{
		return foldedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve smeared plot from unfinalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
}
TH1F * XvsYNormalisedFolding::UncorrectedHistogram()
{
	if ( finalised )
	{
		return inputDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve plot without smearing from unfinalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
}

TH1F * XvsYNormalisedFolding::MCTruthHistogram()
{
	if ( finalised )
	{
		return reconstructedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve MC truth plot from unfinalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
}

//Return a distribution for use in the cross-checks
Distribution * XvsYNormalisedFolding::MonteCarloTruthForCrossCheck() 
{
	//This doesn't really make sense to use
	return XvsYFolder->GetTruthDistribution();
}


TH2F * XvsYNormalisedFolding::SmearingMatrix()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
}

//Convert from an X*Y linearised distribution to an X vs Y plot
//WARNING: this method deletes the argument object
TH1F * XvsYNormalisedFolding::Delinearise( TH1F * LinearisedDistribution )
{
	//Find the target number of bins
	int binNumber = DistributionIndices->GetBinNumber(0);

	//Make a vector of the de-linearised data
	vector<double> delinearisedDistribution( binNumber, 0.0 );
	vector<int> separateIndices;
	vector<double> centralValues, dataCentralValues;
	for ( int binIndex = 0; binIndex < DistributionIndices->GetBinNumber(); binIndex++ )
	{
		//Work out the delinearised bin index and central value
		separateIndices = DistributionIndices->GetNDimensionalIndex(binIndex);
		centralValues = DistributionIndices->GetCentralValues(separateIndices);

		//Increment the bin in the delinearised distribution
		double thisBinValue = LinearisedDistribution->GetBinContent(binIndex) * centralValues[1];
		delinearisedDistribution[ separateIndices[0] ] += thisBinValue;
	}

	//Delete the old distribution
	string name = LinearisedDistribution->GetName();
	string title = LinearisedDistribution->GetTitle();
	delete LinearisedDistribution;

	//Make the new distribution
	TH1F * delinearisedHistogram = new TH1F( name.c_str(), title.c_str(), binNumber - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0] );

	//Copy the data into the new distribution
	for ( int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		delinearisedHistogram->SetBinContent( binIndex, delinearisedDistribution[binIndex] );
	}

	//Return the new distribution
	return delinearisedHistogram;
}

vector<double> XvsYNormalisedFolding::DelineariseErrors( vector<double> LinearisedErrors )
{
	//Find the target number of bins
	int binNumber = DistributionIndices->GetBinNumber(0);

	//Make a vector of the de-linearised data
	vector<double> delinearisedErrors( binNumber, 0.0 );
	vector<int> separateIndices;
	for ( int binIndex = 0; binIndex < DistributionIndices->GetBinNumber(); binIndex++ )
	{
		//Find the delinearised bin index
		separateIndices = DistributionIndices->GetNDimensionalIndex(binIndex);

		//Increment the bin error
		double thisBinValue = LinearisedErrors[binIndex];
		delinearisedErrors[ separateIndices[0] ] += thisBinValue;
	}

	return delinearisedErrors;
}

string XvsYNormalisedFolding::Description( bool WithSpaces )
{
	if (WithSpaces)
	{
		return xName + " vs " + yName;
	}
	else
	{
		return xName + "vs" + yName;
	}
}
string XvsYNormalisedFolding::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector<double> XvsYNormalisedFolding::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedInputErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised XvsYNormalisedFolding" << endl;
		exit(1);
	}
}
