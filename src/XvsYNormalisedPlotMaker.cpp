#include "XvsYNormalisedPlotMaker.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "TFile.h"

using namespace std;

//Default constructor - useless
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker()
{
}

//Constructor with the names to use for the variables
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
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
	XUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + priorName, uniqueID );

	//Store the y range
	minima.push_back( YMinimum );
	maxima.push_back( YMaximum );
	binNumbers.push_back( YBinNumber );

	//Make the x vs y unfolder
	XvsYUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + "vs" + yName + priorName, uniqueID );

	//Set up the indices for the distributions
	DistributionIndices = new DataIndices( binNumbers, minima, maxima );

	//Set up the cross-check for data loss in delinearisation
	char idString[20];
	sprintf( idString, "%d", rand() );
	string xvsyTruthName = xName + yName + priorName + "TruthCheck" + idString;
	string xTruthName = xName + priorName + "TruthCheck" + idString;
	xvsyTruthCheck = new TH1F( xvsyTruthName.c_str(), xvsyTruthName.c_str(), XBinNumber, XMinimum, XMaximum );
	xTruthCheck = new TH1F( xTruthName.c_str(), xTruthName.c_str(), XBinNumber, XMinimum, XMaximum );
}

//Destructor
XvsYNormalisedPlotMaker::~XvsYNormalisedPlotMaker()
{
	delete xvsyTruthCheck;
	delete xTruthCheck;
	delete DistributionIndices;
	delete XUnfolder;
	delete XvsYUnfolder;
	delete correctedDistribution;
	delete uncorrectedDistribution;
	delete mcTruthDistribution;
	delete smearingMatrix;
}

//Copy the object
IPlotMaker * XvsYNormalisedPlotMaker::Clone( string NewPriorName )
{
	return new XvsYNormalisedPlotMaker( xName, yName, NewPriorName,
			DistributionIndices->GetBinNumber(0) - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0],
			DistributionIndices->GetBinNumber(1) - 2, DistributionIndices->GetMinima()[1], DistributionIndices->GetMaxima()[1],
			scaleFactor, uniqueID );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XvsYNormalisedPlotMaker::StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XvsYNormalisedPlotMaker" << endl;
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
		XUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		//Store the y values
		truthValues.push_back( yTruthValue );
		reconstructedValues.push_back( yReconstructedValue );
		XvsYUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		if ( useInPrior )
		{
			//Store values for checking the delinearisation
			xTruthCheck->Fill( xTruthValue, truthWeight );
			xvsyTruthCheck->Fill( xTruthValue, yTruthValue * truthWeight );
		}
	}
}
void XvsYNormalisedPlotMaker::StoreMiss( InputNtuple * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XvsYNormalisedPlotMaker" << endl;
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
		XUnfolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );

		//Store the y value
		truthValues.push_back( yTruthValue );
		XvsYUnfolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );

		if ( useInPrior )
		{
			//Store values for checking the delinearisation
			xTruthCheck->Fill( xTruthValue, truthWeight );
			xvsyTruthCheck->Fill( xTruthValue, yTruthValue * truthWeight );
		}
	}
}
void XvsYNormalisedPlotMaker::StoreFake( InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XvsYNormalisedPlotMaker" << endl;
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
		XUnfolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );

		//Store the y value
		reconstructedValues.push_back( yReconstructedValue );
		XvsYUnfolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );
	}
}
void XvsYNormalisedPlotMaker::StoreData( InputNtuple * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XvsYNormalisedPlotMaker" << endl;
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
		XUnfolder->StoreDataValue( dataValues, dataWeight );

		//Store the y value
		dataValues.push_back( yDataValue );
		XvsYUnfolder->StoreDataValue( dataValues, dataWeight );

		//Store values for performing the delinearisation
		DistributionIndices->StoreDataValue( dataValues, dataWeight );
	} 
}

//Do the unfolding
void XvsYNormalisedPlotMaker::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XvsYNormalisedPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Closure tests
		XUnfolder->ClosureTest();
		XvsYUnfolder->ClosureTest();

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
		TH1F * DelinearisedXvsYCorrected = Delinearise(XvsYCorrected);
		TH1F * DelinearisedXvsYUncorrected = Delinearise(XvsYUncorrected);
		TH1F * DelinearisedXvsYTruth = Delinearise(XvsYTruth);

		//Get the error vectors
		vector<double> XErrors = XUnfolder->SumOfDataWeightSquares();
		vector<double> XvsYErrors = XvsYUnfolder->SumOfDataWeightSquares();

		//Delinearise the x vs y errors
		vector<double> delinearisedXvsYErrors = DelineariseErrors( XvsYErrors );

		//Combine errors
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			//Add the errors from the x vs y distribution and the divisor
			//The formula is stright from ROOT::TH1::Divide, so I hope it's right
			double XBinValue = XCorrected->GetBinContent(binIndex);
			double XvsYBinValue = DelinearisedXvsYCorrected->GetBinContent(binIndex);
			double componentOne = XErrors[binIndex] * XvsYBinValue * XvsYBinValue;
			double componentTwo = delinearisedXvsYErrors[binIndex] * XBinValue * XBinValue;
			double componentThree = XBinValue * XBinValue;
			double combinedError = ( sqrt( componentOne + componentTwo ) / componentThree );
			combinedError *= scaleFactor;
			correctedDataErrors.push_back( combinedError );
		}

		//Normalise the x vs y distributions and scale appropriately
		DelinearisedXvsYCorrected->Divide( DelinearisedXvsYCorrected, XCorrected, scaleFactor, 1.0 );
		DelinearisedXvsYUncorrected->Divide( DelinearisedXvsYUncorrected, XUncorrected, scaleFactor, 1.0 );
		DelinearisedXvsYTruth->Divide( DelinearisedXvsYTruth, XTruth, scaleFactor, 1.0 );

		//Check for data loss in the delinearisation
		xvsyTruthCheck->Divide( xvsyTruthCheck, xTruthCheck, scaleFactor, 1.0 );
		bool dataLost = false;
		double averagePercentError = 0.0;
		for ( int binIndex = 0; binIndex <= xTruthCheck->GetNbinsX() + 1; binIndex++ )
		{
			//Compare the delinearised value with one calculated without going through delinearisation
			double correctValue = xvsyTruthCheck->GetBinContent( binIndex );
			double delinearisedValue = DelinearisedXvsYTruth->GetBinContent( binIndex );
			double percentError = ( delinearisedValue - correctValue ) * 100.0 / correctValue;

			//Check for stupid values
			if ( isnan( percentError ) )
			{
				percentError = 0.0;
			}

			//Increment average error
			averagePercentError += fabs( percentError );

			//If there's a > 1%  discrepancy between the delinearised truth value and the correct value, warn that data is being lost in binning
			if ( fabs( percentError ) > 1.0 )
			{
				cerr << "Bin " << binIndex << " has value " << delinearisedValue << " which has a " << percentError << "\% error vs the reference value (" << correctValue << ")" << endl;
				dataLost = true;
			}
		}
		averagePercentError /= (double)( xTruthCheck->GetNbinsX() + 2 );
		cout << "Average bin error from delinearisation: " << averagePercentError << "\%" << endl;
		if ( dataLost || averagePercentError > 0.5 )
		{
			cerr << "Delinearisation has caused significant changes in distribution bins compared to their reference values" << endl;
			cerr << "Improve the binning of " << yName << " in the " << xName << " vs " << yName << " plot" << endl;
			cerr << "Try using a finer binning, and avoid large numbers of events in under or overflow bins" << endl;
			exit(1);
		}

		//Free some memory
		delete XCorrected;
		delete XUncorrected;
		delete XTruth;

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

//Make a cross-check with MC
int XvsYNormalisedPlotMaker::MonteCarloCrossCheck( TH1F * ReferencePlot, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	return XvsYUnfolder->MonteCarloCrossCheck( ReferencePlot, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
	//return XUnfolder->MonteCarloCrossCheck( ReferencePlot, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
}

//Return some plots
TH1F * XvsYNormalisedPlotMaker::CorrectedDistribution()
{
	if ( finalised )
	{
		return correctedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve corrected plot from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XvsYNormalisedPlotMaker::UncorrectedDistribution()
{
	if ( finalised )
	{
		return uncorrectedDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve uncorrected plot from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XvsYNormalisedPlotMaker::MCTruthDistribution()
{
	if ( finalised )
	{
		return mcTruthDistribution;
	}       
	else
	{
		//Make some plot titles
		char uniqueIDString[10];
		sprintf( uniqueIDString, "%d", uniqueID );
		string XvsYName = xName + "vs" + yName + priorName + uniqueIDString;
		string XvsYTitle = xName + " vs " + yName + " using " + priorName;

		//Retrieve the results
		return XvsYUnfolder->GetTruthDistribution( XvsYName + "RawTruth", XvsYTitle + " Raw Truth Distribution" );
		//return XUnfolder->GetTruthDistribution( XvsYName + "RawTruth", XvsYTitle + " Raw Truth Distribution" );
	}
}
TH2F * XvsYNormalisedPlotMaker::SmearingMatrix()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}

//Convert from an X*Y linearised distribution to an X vs Y plot
//WARNING: this method deletes the argument object
TH1F * XvsYNormalisedPlotMaker::Delinearise( TH1F * LinearisedDistribution )
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

vector<double> XvsYNormalisedPlotMaker::DelineariseErrors( vector<double> LinearisedErrors )
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

string XvsYNormalisedPlotMaker::Description( bool WithSpaces )
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
string XvsYNormalisedPlotMaker::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector<double> XvsYNormalisedPlotMaker::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedDataErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}

