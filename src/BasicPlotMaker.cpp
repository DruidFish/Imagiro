#include "BasicPlotMaker.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "TFile.h"

using namespace std;

//Default constructor - useless
BasicPlotMaker::BasicPlotMaker()
{
}

//Constructor with the names to use for the variables
BasicPlotMaker::BasicPlotMaker( string XVariableName, string YVariableName, string PriorName,
		int XBinNumber, double XMinimum, double XMaximum,
		double ScaleFactor, int UniqueID ) : xName( XVariableName ), yName( YVariableName ), priorName( PriorName ), finalised( false ), uniqueID(UniqueID), scaleFactor(ScaleFactor)
{
	vector<double> minima, maxima;
	vector<int> binNumbers;

	//Store the x range
	minima.push_back( XMinimum );
	maxima.push_back( XMaximum );
	binNumbers.push_back( XBinNumber );

	//Make the x unfolder
	XUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + priorName, uniqueID );

	//Make the x vs y unfolder
	XvsYUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + "vs" + yName + priorName, uniqueID );

	//Set up the indices for the distributions
	DistributionIndices = new Indices( binNumbers, minima, maxima );
}

//Destructor
BasicPlotMaker::~BasicPlotMaker()
{
	delete DistributionIndices;
	delete XUnfolder;
	delete XvsYUnfolder;
	delete correctedDistribution;
	delete uncorrectedDistribution;
	delete mcTruthDistribution;
	delete smearingMatrix;
}

//Copy the object
IPlotMaker * BasicPlotMaker::Clone( string NewPriorName )
{
	return new BasicPlotMaker( xName, yName, NewPriorName,
			DistributionIndices->GetBinNumber(0) - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0],
			scaleFactor, uniqueID );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void BasicPlotMaker::StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised BasicPlotMaker" << endl;
		exit(1);
	}
	else
	{
		vector<double> truthValues, reconstructedValues;

		//Get the combined weight
		double truthWeight = TruthInput->EventWeight();
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Store the x values
		truthValues.push_back( TruthInput->GetValue( xName ) );
		reconstructedValues.push_back( ReconstructedInput->GetValue( xName ) );
		XUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		//Store the y values
		truthWeight *= TruthInput->GetValue( yName );
		XvsYUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );
	}
}
void BasicPlotMaker::StoreMiss( InputNtuple * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised BasicPlotMaker" << endl;
		exit(1);
	}
	else
	{
		vector<double> truthValues;

		//Get the weight
		double weight = TruthInput->EventWeight();

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Store the x value
		truthValues.push_back( TruthInput->GetValue( xName ) );
		XUnfolder->StoreUnreconstructedTruth( truthValues, weight, useInPrior );

		//Store the y value
		weight *= TruthInput->GetValue( yName );
		XvsYUnfolder->StoreUnreconstructedTruth( truthValues, weight, useInPrior );
	}
}
void BasicPlotMaker::StoreFake( InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised BasicPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector<double> reconstructedValues;

		//Get the weight
		double weight = ReconstructedInput->EventWeight();

		//Store the x value
		reconstructedValues.push_back( ReconstructedInput->GetValue( xName ) );
		XUnfolder->StoreReconstructedFake( reconstructedValues, weight );

		//Store the y value
		weight *= ReconstructedInput->GetValue( yName );
		XvsYUnfolder->StoreReconstructedFake( reconstructedValues, weight );
	}
}
void BasicPlotMaker::StoreData( InputNtuple * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised BasicPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector<double> dataValues;

		//Get the weight
		double weight = DataInput->EventWeight();

		//Store the x value
		dataValues.push_back( DataInput->GetValue( xName ) );
		XUnfolder->StoreDataValue( dataValues, weight );

		//Store the y value
		weight *= DataInput->GetValue( yName );
		XvsYUnfolder->StoreDataValue( dataValues, weight );
	} 
}

//Do the unfolding
void BasicPlotMaker::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "BasicPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Unfold the distributions
		XUnfolder->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
		XvsYUnfolder->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );

		//Make some plot titles
		char uniqueIDString[10];
		sprintf( uniqueIDString, "%d", uniqueID );
		string XFullName = xName + priorName + uniqueIDString;
		string XFullTitle = xName + " using " + priorName;
		string XvsYName = xName + "vs" + yName + priorName + uniqueIDString;
		string XvsYTitle = xName + " vs " + yName + " using " + priorName;

		//Retrieve the results
		TH1F * XCorrected = XUnfolder->UnfoldedDistribution( XFullName + "Corrected", XFullTitle + " Corrected Distribution" );
		TH1F * XvsYCorrected = XvsYUnfolder->UnfoldedDistribution( XvsYName + "Corrected", XvsYTitle + " Corrected Distribution" );

		//Retrieve some other bits for debug
		TH1F * XUncorrected = XUnfolder->GetUncorrectedDataDistribution( XFullName + "Uncorrected", XFullTitle + " Uncorrected Distribution" );
		TH1F * XvsYUncorrected = XvsYUnfolder->GetUncorrectedDataDistribution( XvsYName + "Uncorrected", XvsYTitle + " Uncorrected Distribution" );

		TH1F * XTruth = XUnfolder->GetTruthDistribution( XFullName + "Truth", XFullTitle + " Truth Distribution" );
		TH1F * XvsYTruth = XvsYUnfolder->GetTruthDistribution( XvsYName + "Truth", XvsYTitle + " Truth Distribution" );

		//TH2F * XSmearing = XUnfolder->GetSmearingMatrix( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );
		TH2F * XvsYSmearing = XvsYUnfolder->GetSmearingMatrix( XvsYName + "Smearing", XvsYTitle + " Smearing Matrix" );

		//De-linearise the x vs y distributions
		TH1F * DelinearisedXvsYCorrected = XvsYCorrected;
		TH1F * DelinearisedXvsYUncorrected = XvsYUncorrected;
		TH1F * DelinearisedXvsYTruth = XvsYTruth;

		//Get the error vectors
		vector<double> XErrors = XUnfolder->SumOfDataWeightSquares();
		vector<double> XvsYErrors = XvsYUnfolder->SumOfDataWeightSquares();

		//Combine errors
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			//Add the errors from the x vs y distribution and the divisor
			//The formula is stright from ROOT::TH1::Divide, so I hope it's right
			double XBinValue = XCorrected->GetBinContent(binIndex);
			double XvsYBinValue = DelinearisedXvsYCorrected->GetBinContent(binIndex);
			double componentOne = XErrors[binIndex] * XvsYBinValue * XvsYBinValue;
			double componentTwo = XvsYErrors[binIndex] * XBinValue * XBinValue;
			double componentThree = XBinValue * XBinValue;
			double combinedError = ( sqrt( componentOne + componentTwo ) / componentThree );
			combinedError *= scaleFactor;
			correctedDataErrors.push_back( combinedError );
		}

		//Normalise the x vs y distributions and scale appropriately
		DelinearisedXvsYCorrected->Divide( DelinearisedXvsYCorrected, XCorrected, scaleFactor, 1.0 );
		DelinearisedXvsYUncorrected->Divide( DelinearisedXvsYUncorrected, XUncorrected, scaleFactor, 1.0 );
		DelinearisedXvsYTruth->Divide( DelinearisedXvsYTruth, XTruth, scaleFactor, 1.0 );

		//Get the y range to plot
		double yMinimum = DistributionIndices->GetMinima()[1] * scaleFactor;
		double yMaximum = DistributionIndices->GetMaxima()[1] * scaleFactor;

		//Format and save the corrected distribution
		correctedDistribution = DelinearisedXvsYCorrected;
		correctedDistribution->SetXTitle( xName.c_str() );
		correctedDistribution->SetYTitle( yName.c_str() );
		correctedDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the uncorrected distribution
		uncorrectedDistribution = DelinearisedXvsYUncorrected;
		uncorrectedDistribution->SetXTitle( xName.c_str() );
		uncorrectedDistribution->SetYTitle( yName.c_str() );
		uncorrectedDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the truth distribution
		mcTruthDistribution = DelinearisedXvsYTruth;
		mcTruthDistribution->SetXTitle( xName.c_str() );
		mcTruthDistribution->SetYTitle( yName.c_str() );
		mcTruthDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the smearing matrix
		smearingMatrix = XvsYSmearing;
		string smearingXTitle = xName + " vs " + yName + " Truth";
		string smearingYTitle = xName + " vs " + yName + " Reconstructed";
		smearingMatrix->SetXTitle( smearingXTitle.c_str() );
		smearingMatrix->SetYTitle( smearingYTitle.c_str() );

		//Bin-by-bin scaling of errors using the corrected data
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double errorScaleFactor = correctedDistribution->GetBinContent(binIndex) / uncorrectedDistribution->GetBinContent(binIndex);
			correctedDataErrors[binIndex] *= errorScaleFactor;
		}

		//Mark as done
		finalised = true;
	}
}

//Return some plots
TH1F * BasicPlotMaker::CorrectedDistribution()
{
	if ( finalised )
	{
		return correctedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve corrected plot from unfinalised BasicPlotMaker" << endl;
		exit(1);
	}
}
TH1F * BasicPlotMaker::UncorrectedDistribution()
{
	if ( finalised )
	{
		return uncorrectedDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve uncorrected plot from unfinalised BasicPlotMaker" << endl;
		exit(1);
	}
}
TH1F * BasicPlotMaker::MCTruthDistribution()
{
	if ( finalised )
	{
		return mcTruthDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve mc truth plot from unfinalised BasicPlotMaker" << endl;
		exit(1);
	}
}
TH2F * BasicPlotMaker::SmearingMatrix()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised BasicPlotMaker" << endl;
		exit(1);
	}
}

string BasicPlotMaker::Description( bool WithSpaces )
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
string BasicPlotMaker::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector<double> BasicPlotMaker::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedDataErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised BasicPlotMaker" << endl;
		exit(1);
	}
}

