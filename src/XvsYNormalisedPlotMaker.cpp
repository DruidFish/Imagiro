/**
  @class XvsYNormalisedPlotMaker

  Unfolds a 2D distribution, and divides it by the unfolded 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#include "XvsYNormalisedPlotMaker.h"
#include "BayesianUnfolding.h"
#include "Folding.h"
#include "NoCorrection.h"
#include "BinByBinUnfolding.h"
#include "UniformIndices.h"
#include "CustomIndices.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

const int BAYESIAN_MODE = 2;

//Default constructor - useless
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker()
{
}

//Constructor with the names to use for the variables
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
		unsigned int XBinNumber, double XMinimum, double XMaximum,
		unsigned int YBinNumber, double YMinimum, double YMaximum, int CorrectionMode, double ScaleFactor )
{
	correctionType = CorrectionMode;
	xName = XVariableName;
	yName = YVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	vector< double > minima, maxima;
	vector< unsigned int > binNumbers;
	systematicRandom = 0;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Store the x range
	minima.push_back( XMinimum );
	maxima.push_back( XMaximum );
	binNumbers.push_back( XBinNumber );

	//Store the y range
	minima.push_back( YMinimum );
	maxima.push_back( YMaximum );
	binNumbers.push_back( YBinNumber );

	//Make the x vs y unfolder
	distributionIndices = new UniformIndices( binNumbers, minima, maxima );
	XvsYUnfolder = MakeCorrector( correctionType, distributionIndices, xName + "vs" + yName + priorName, thisPlotID );

	//Set up the cross-check for data loss in delinearisation
	stringstream idString;
	idString << thisPlotID;
	string xvsyTruthName = xName + yName + priorName + "TruthCheck" + idString.str();
	xvsyTruthCheck = new TProfile( xvsyTruthName.c_str(), xvsyTruthName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );

	//Make a summary for the y data values
	yValueSummary = new StatisticsSummary();

	//Avoid killing ROOT
	doPlotSmearing = ( XBinNumber * YBinNumber < 1200 );

	//Error calculation
	xvsyTruthName += "SimpleProfile";
	simpleDataProfile = new TProfile( xvsyTruthName.c_str(), xvsyTruthName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );
	simpleDataProfile->Sumw2();
}

//Constructor with the names to use for the variables
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
		vector< double > XBinLowEdges, vector< double > YBinLowEdges, int CorrectionMode, double ScaleFactor )
{
	correctionType = CorrectionMode;
	xName = XVariableName;
	yName = YVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	vector< vector< double > > binEdges;
	systematicRandom = 0;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Store the x range
	binEdges.push_back( XBinLowEdges );

	//Store the y range
	binEdges.push_back( YBinLowEdges );

	//Make the x vs y unfolder
	distributionIndices = new CustomIndices( binEdges );
	XvsYUnfolder = MakeCorrector( correctionType, distributionIndices, xName + "vs" + yName + priorName, thisPlotID );

	//Set up the cross-check for data loss in delinearisation
	stringstream idString;
	idString << thisPlotID;
	string xvsyTruthName = xName + yName + priorName + "TruthCheck" + idString.str();
	xvsyTruthCheck = new TProfile( xvsyTruthName.c_str(), xvsyTruthName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );

	//Make a summary for the y data values
	yValueSummary = new StatisticsSummary();

	//Avoid killing ROOT
	doPlotSmearing = ( distributionIndices->GetBinNumber(0) * distributionIndices->GetBinNumber(1) < 1200 );

	//Error calculation
	xvsyTruthName += "SimpleProfile";
	simpleDataProfile = new TProfile( xvsyTruthName.c_str(), xvsyTruthName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );
	simpleDataProfile->Sumw2();
}

//To be used only with Clone
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
		IIndexCalculator * DistributionIndices, int CorrectionMode, unsigned int OriginalID, double ScaleFactor, vector< vector< double > > InputOffsets, vector< vector< double > > InputWidths )
{
	correctionType = CorrectionMode;
	xName = XVariableName;
	yName = YVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;

	systematicOffsets = InputOffsets;
	systematicWidths = InputWidths;
	if ( systematicWidths.size() == 0 )
	{
		systematicRandom = 0;
	}
	else
	{
		systematicRandom = new TRandom3(0);
	}

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID + OriginalID;

	//Make the x vs y unfolder
	distributionIndices = DistributionIndices;
	XvsYUnfolder = MakeCorrector( correctionType, distributionIndices, xName + "vs" + yName + priorName, thisPlotID );

	//Make the systematics too
	for ( unsigned int experimentIndex = 0; experimentIndex < systematicWidths.size(); experimentIndex++ )
	{
		systematicUnfolders.push_back( XvsYUnfolder->CloneShareSmearingMatrix() );
	}

	//Set up the cross-check for data loss in delinearisation
	stringstream idString;
	idString << thisPlotID;
	string xvsyTruthName = xName + yName + priorName + "TruthCheck" + idString.str();
	xvsyTruthCheck = new TProfile( xvsyTruthName.c_str(), xvsyTruthName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );

	//Make a summary for the y data values
	yValueSummary = new StatisticsSummary();

	//Avoid killing ROOT
	doPlotSmearing = ( distributionIndices->GetBinNumber(0) * distributionIndices->GetBinNumber(1) < 1200 );

	//Error calculation
	xvsyTruthName += "SimpleProfile";
	simpleDataProfile = new TProfile( xvsyTruthName.c_str(), xvsyTruthName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );
	simpleDataProfile->Sumw2();
}

//Destructor
XvsYNormalisedPlotMaker::~XvsYNormalisedPlotMaker()
{
	delete xvsyTruthCheck;
	delete distributionIndices;
	delete XvsYUnfolder;
	for ( unsigned int experimentIndex = 0; experimentIndex < systematicUnfolders.size(); experimentIndex++ )
	{
		delete systematicUnfolders[ experimentIndex ];
	}
	delete simpleDataProfile;
	if ( finalised )
	{
		delete correctedDistribution;
		delete uncorrectedDistribution;
		delete mcTruthDistribution;
		if ( smearingMatrix )
		{
			delete smearingMatrix;
		}
		if ( covarianceMatrix )
		{
			delete covarianceMatrix;
		}

		for ( unsigned int experimentIndex = 0; experimentIndex < systematicResults.size(); experimentIndex++ )
		{
			delete systematicResults[ experimentIndex ];
		}
	}
	if ( systematicRandom )
	{
		delete systematicRandom;
	}
}

//Copy the object
XvsYNormalisedPlotMaker * XvsYNormalisedPlotMaker::Clone( string NewPriorName )
{
	return new XvsYNormalisedPlotMaker( xName, yName, NewPriorName, distributionIndices->Clone(), correctionType, thisPlotID, scaleFactor, systematicOffsets, systematicWidths );
}

//Set up a systematic error study
void XvsYNormalisedPlotMaker::AddSystematic( vector< double > SystematicOffset, vector< double > SystematicWidth, unsigned int NumberOfPseudoExperiments )
{
	//Check for correct number of observables
	if ( SystematicOffset.size() != 2 || SystematicWidth.size() != 2 )
	{
		cerr << "ERROR: There's two observables in an XvsYNormalisedPlotmaker, please define systematic errors on both, even if one is just 0.0"<< endl;
		exit(1);
	}

	//Check for sensible pseudoexperiment number
	double experimentNumber = NumberOfPseudoExperiments;
	if ( SystematicWidth[0] == 0.0 && SystematicWidth[1] == 0.0 && NumberOfPseudoExperiments > 1 )
	{
		cout << "INFO: There's no point doing many pseudoexperiments for just a systematic offset - using 1" << endl;
		experimentNumber = 1;
	}
	if ( ( SystematicWidth[0] != 0.0 || SystematicWidth[1] != 0.0 ) && NumberOfPseudoExperiments < 1 )
	{
		cerr << "ERROR: You need to pick a number of pseudoexperiments to propagate this systematic width" << endl;
		exit(1);
	}

	//Check for dodgy business
	if ( ( SystematicWidth[0] != 0.0 && SystematicOffset[0] != 0.0 ) || ( SystematicWidth[1] != 0.0 && SystematicOffset[1] != 0.0 ) )
	{
		cout << "WARNING: Doing a systematic width and offset at the same time can bias the final result. Tread carefully!" << endl;
	}

	//Make the random number generator
	if ( !systematicRandom && ( SystematicWidth[0] != 0.0 || SystematicWidth[1] != 0.0 ) )
	{
		systematicRandom = new TRandom3(0);
	}

	//Make the extra experiments
	for ( unsigned int experimentIndex = 0; experimentIndex < experimentNumber; experimentIndex++ )
	{
		systematicUnfolders.push_back( XvsYUnfolder->CloneShareSmearingMatrix() );
		systematicOffsets.push_back( SystematicOffset );
		systematicWidths.push_back( SystematicWidth );
	}
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XvsYNormalisedPlotMaker::StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
	else
	{
		vector< double > truthValues, reconstructedValues;

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

		//Store the y values
		truthValues.push_back( yTruthValue );
		reconstructedValues.push_back( yReconstructedValue );
		XvsYUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		if ( useInPrior )
		{
			//Store values for checking the delinearisation
			xvsyTruthCheck->Fill( xTruthValue, yTruthValue, truthWeight );
		}
	}
}
void XvsYNormalisedPlotMaker::StoreMiss( IFileInput * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
	else
	{
		vector< double > truthValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xTruthValue = TruthInput->GetValue( xName );
		double yTruthValue = TruthInput->GetValue( yName );
		double truthWeight = TruthInput->EventWeight();

		//Store the x value
		truthValues.push_back( xTruthValue );

		//Store the y value
		truthValues.push_back( yTruthValue );
		XvsYUnfolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );

		if ( useInPrior )
		{
			//Store values for checking the delinearisation
			xvsyTruthCheck->Fill( xTruthValue, yTruthValue, truthWeight );
		}
	}
}
void XvsYNormalisedPlotMaker::StoreFake( IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector< double > reconstructedValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( ReconstructedInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xReconstructedValue = ReconstructedInput->GetValue( xName );
		double yReconstructedValue = ReconstructedInput->GetValue( yName );
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x value
		reconstructedValues.push_back( xReconstructedValue );

		//Store the y value
		reconstructedValues.push_back( yReconstructedValue );
		XvsYUnfolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );
	}
}
void XvsYNormalisedPlotMaker::StoreData( IFileInput * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		double randomOne, randomTwo;
		vector<double> dataValues;

		//Retrieve the values from the Ntuple
		double xDataValue = DataInput->GetValue( xName );
		double yDataValue = DataInput->GetValue( yName );
		double dataWeight = DataInput->EventWeight();

		//Store the x value
		dataValues.push_back( xDataValue );

		//Store the y value
		dataValues.push_back( yDataValue );
		XvsYUnfolder->StoreDataValue( dataValues, dataWeight );
		yValueSummary->StoreEvent( yDataValue, dataWeight );

		//Store values for performing the delinearisation
		distributionIndices->StoreDataValue( dataValues, dataWeight );

		//Store the values for error calculation
		simpleDataProfile->Fill( xDataValue, yDataValue, dataWeight );

		//Do the systematic experiments too
		for ( unsigned int experimentIndex = 0; experimentIndex < systematicOffsets.size(); experimentIndex++ )
		{
			dataValues[0] = xDataValue + systematicOffsets[ experimentIndex ][0];
			dataValues[1] = yDataValue + systematicOffsets[ experimentIndex ][1];

			if ( systematicWidths[ experimentIndex ][0] != 0.0 || systematicWidths[ experimentIndex ][1] != 0.0 )
			{
				systematicRandom->Rannor( randomOne, randomTwo );
				dataValues[0] += ( randomOne * systematicWidths[ experimentIndex ][0] );
				dataValues[1] += ( randomTwo * systematicWidths[ experimentIndex ][1] );
			}

			systematicUnfolders[ experimentIndex ]->StoreDataValue( dataValues, dataWeight );
		}
	}
}

//Do the unfolding
void XvsYNormalisedPlotMaker::Correct( unsigned int MostIterations, bool SkipUnfolding, unsigned int ErrorMode, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XvsYNormalisedPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Unfold the distributions
		if ( !SkipUnfolding )
		{
			XvsYUnfolder->Correct( MostIterations, ErrorMode, WithSmoothing );

			//Systematics
			for ( unsigned int experimentIndex = 0; experimentIndex < systematicUnfolders.size(); experimentIndex++ )
			{
				systematicUnfolders[ experimentIndex ]->Correct( MostIterations, ErrorMode, WithSmoothing );
			}
		}

		//Make some plot titles
		stringstream uniqueIDString;
		uniqueIDString << thisPlotID;
		string XFullName = xName + priorName + uniqueIDString.str();
		string XFullTitle = xName + " using " + priorName;
		string XvsYName = xName + "vs" + yName + priorName + uniqueIDString.str();
		string XvsYTitle = xName + " vs " + yName + " using " + priorName;

		//Retrieve the results
		TH1F * XvsYCorrected = XvsYUnfolder->GetCorrectedHistogram( XvsYName + "Corrected", XvsYTitle + " Corrected Distribution" );

		//Retrieve some other bits for debug
		TH1F * XvsYUncorrected = XvsYUnfolder->GetUncorrectedHistogram( XvsYName + "Uncorrected", XvsYTitle + " Uncorrected Distribution" );
		TH1F * XvsYTruth = XvsYUnfolder->GetTruthHistogram( XvsYName + "Truth", XvsYTitle + " Truth Distribution" );
		TH1F * XvsYReco = XvsYUnfolder->GetReconstructedHistogram( XvsYName + "Reco", XvsYTitle + " Reco Distribution" );

		//De-linearise the x vs y distributions
		TH1F * DelinearisedXvsYCorrected = MakeProfile( XvsYCorrected );
		TH1F * DelinearisedXvsYUncorrected = MakeProfile( XvsYUncorrected );
		TH1F * DelinearisedXvsYTruth = MakeProfile( XvsYTruth );
		TH1F * DelinearisedXvsYReco = MakeProfile( XvsYReco );

		//And the systematics
		for ( unsigned int experimentIndex = 0; experimentIndex < systematicUnfolders.size(); experimentIndex++ )
		{
			//Name
			stringstream systematicString;
			systematicString << thisPlotID << XvsYName << priorName << "Systematic" << experimentIndex;

			//Rerieve
			TH1F * systematicPlot = systematicUnfolders[ experimentIndex ]->GetCorrectedHistogram( systematicString.str(), systematicString.str() );

			//Delinearise
			systematicPlot = MakeProfile( systematicPlot );

			//Scale
			systematicPlot->Scale( scaleFactor );

			//Store
			systematicResults.push_back( systematicPlot );
		}

		//Make a vector of bin error values
		if ( correctionType != BAYESIAN_MODE || ErrorMode < 1 )
		{
			for ( unsigned int binIndex = 0; binIndex < distributionIndices->GetBinNumber(0); binIndex++ )
			{
				//Get the error calculated by the uncorrected data TProfile, then scale it just like the bin values
				double combinedError = simpleDataProfile->GetBinError( binIndex ) * scaleFactor;

				//Store error
				correctedDataErrors.push_back( combinedError );
			}
		}
		else
		{
			//Get the errors from the unfolder, then delinearise
			vector<double> XvsYErrors = XvsYUnfolder->Variances();
			vector<double> delinearisedXvsYErrors = DelineariseErrors( XvsYErrors );

			for ( unsigned int binIndex = 0; binIndex < distributionIndices->GetBinNumber(0); binIndex++ )
			{
				//Scale the error, just like the bin value
				double combinedError = delinearisedXvsYErrors[ binIndex ] * scaleFactor;

				//Store error
				correctedDataErrors.push_back( combinedError );
			}
		}

		//Scale the bin values
		DelinearisedXvsYCorrected->Scale( scaleFactor );
		DelinearisedXvsYUncorrected->Scale( scaleFactor );
		DelinearisedXvsYTruth->Scale( scaleFactor );
		DelinearisedXvsYReco->Scale( scaleFactor );
		xvsyTruthCheck->Scale( scaleFactor );

		//Check for data loss in the delinearisation
		bool dataLost = false;
		double averagePercentError = 0.0;
		for ( int binIndex = 0; binIndex < xvsyTruthCheck->GetNbinsX() + 2; binIndex++ )
		{
			//Compare the delinearised value with one calculated without going through delinearisation
			double correctValue = xvsyTruthCheck->GetBinContent( binIndex );
			double delinearisedValue = DelinearisedXvsYTruth->GetBinContent( binIndex );
			if ( correctValue != 0.0 )
			{
				double errorFraction = fabs( delinearisedValue - correctValue ) / correctValue;
				double percentError = errorFraction * 100.0;

				//Check for stupid values
				if ( isnan( percentError ) )
				{
					percentError = 0.0;
				}

				//Add this error to the overall error calculation
				correctedDataErrors[ binIndex ] *= ( 1.0 + errorFraction );

				//Increment average error
				averagePercentError += percentError;

				//If there's a > 1%  discrepancy between the delinearised truth value and the correct value, warn that data is being lost in binning
				if ( fabs( percentError ) > 1.0 )
				{
					cerr << "Bin " << binIndex << " has value " << delinearisedValue << " which has a " << percentError << "\% error vs the reference value (" << correctValue << ")" << endl;
					dataLost = true;
				}
			}
		}
		averagePercentError /= (double)( xvsyTruthCheck->GetNbinsX() + 2 );
		cout << "Average bin error from delinearisation: " << averagePercentError << "\%" << endl;
		if ( dataLost || averagePercentError > 0.5 )
		{
			cerr << "Delinearisation has caused significant changes in distribution bins compared to their reference values" << endl;
			cerr << "Improve the binning of " << yName << " in the " << xName << " vs " << yName << " plot" << endl;
			cerr << "Try using a finer binning, and avoid large numbers of events in under or overflow bins" << endl;
			cerr << "Suggested bin width: " << yValueSummary->OptimumBinWidth() << endl;
			//exit(1);
		}

		//Free some memory
		delete yValueSummary;

		//Format and save the corrected distribution
		correctedDistribution = DelinearisedXvsYCorrected;
		correctedDistribution->SetXTitle( xName.c_str() );
		correctedDistribution->SetYTitle( yName.c_str() );

		//Format and save the uncorrected distribution
		uncorrectedDistribution = DelinearisedXvsYUncorrected;
		uncorrectedDistribution->SetXTitle( xName.c_str() );
		uncorrectedDistribution->SetYTitle( yName.c_str() );

		//Format and save the truth distribution
		mcTruthDistribution = DelinearisedXvsYTruth;
		mcTruthDistribution->SetXTitle( xName.c_str() );
		mcTruthDistribution->SetYTitle( yName.c_str() );

		//Format and save the reco distribution
		mcRecoDistribution = DelinearisedXvsYReco;
		mcRecoDistribution->SetXTitle( xName.c_str() );
		mcRecoDistribution->SetYTitle( yName.c_str() );

		//Format and save the smearing matrix
		if ( doPlotSmearing )
		{
			smearingMatrix = XvsYUnfolder->GetSmearingHistogram( XvsYName + "Smearing", XvsYTitle + " Smearing Matrix" );
			string smearingXTitle = xName + " vs " + yName + " Truth";
			string smearingYTitle = xName + " vs " + yName + " Reconstructed";
			smearingMatrix->SetXTitle( smearingXTitle.c_str() );
			smearingMatrix->SetYTitle( smearingYTitle.c_str() );
		}
		else
		{
			smearingMatrix = NULL;
		}

		//Format and save the covariance matrix
		if ( ErrorMode > 1 )
		{
			covarianceMatrix = XvsYUnfolder->DAgostiniCovariance( XvsYName + "Covariance", XvsYTitle + " Covariance Matrix" );
			covarianceMatrix->SetXTitle( "Bin Number" );
			covarianceMatrix->SetYTitle( "Bin Number" );
		}
		else
		{
			covarianceMatrix = NULL;
		}

		//Bin-by-bin scaling of errors using the corrected data if a full error calculation was not done
		if ( correctionType != BAYESIAN_MODE || ErrorMode < 1 )
		{
			for ( unsigned int binIndex = 0; binIndex < correctedDataErrors.size(); binIndex++ )
			{
				double errorScaleFactor = correctedDistribution->GetBinContent( binIndex ) / uncorrectedDistribution->GetBinContent( binIndex );

				//Check for div0 errors
				if ( isinf( errorScaleFactor ) )
				{
					errorScaleFactor = 1.0;
				}
				if ( isnan( errorScaleFactor ) )
				{
					errorScaleFactor = 0.0;
				}

				correctedDataErrors[ binIndex ] *= errorScaleFactor;
			}
		}

		//Mark as done
		finalised = true;
	}
}

//Do a closure test
bool XvsYNormalisedPlotMaker::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
{
	return XvsYUnfolder->ClosureTest( MostIterations, WithSmoothing );
}

//Make a cross-check with MC
unsigned int XvsYNormalisedPlotMaker::MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing )
{
	return XvsYUnfolder->MonteCarloCrossCheck( InputPriorDistribution, InputSmearing, WithSmoothing );
}

//Return some plots
TH1F * XvsYNormalisedPlotMaker::CorrectedHistogram()
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
TH1F * XvsYNormalisedPlotMaker::UncorrectedHistogram()
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
TH1F * XvsYNormalisedPlotMaker::MCTruthHistogram()
{
	if ( finalised )
	{
		return mcTruthDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve MC truth plot from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XvsYNormalisedPlotMaker::MCRecoHistogram()
{
	if ( finalised )
	{
		return mcRecoDistribution;
	}
	else
	{
		cerr << "Trying to retrieve MC reco plot from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}


//Return a distribution for use in the cross-checks
Distribution * XvsYNormalisedPlotMaker::PriorDistributionForCrossCheck()
{
	return XvsYUnfolder->GetTruthDistribution();
}
SmearingMatrix * XvsYNormalisedPlotMaker::SmearingMatrixForCrossCheck()
{
	return XvsYUnfolder->GetSmearingMatrix();
}

TH2F * XvsYNormalisedPlotMaker::SmearingHistogram()
{
	if ( finalised )
	{
		if ( doPlotSmearing )
		{
			return smearingMatrix;
		}
		else
		{
			cout << "Not plotting smearing matrix: it's too big, and Root won't like it" << endl;
			return NULL;
		}
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}

//Return the results of the systematic experiments
vector< TH1F* > XvsYNormalisedPlotMaker::SystematicHistograms()
{
	return systematicResults;
}

//Convert from an X*Y linearised distribution to an X vs Y plot
//WARNING: this method deletes the argument object
TH1F * XvsYNormalisedPlotMaker::MakeProfile( TH1F * LinearisedDistribution )
{
	//Make an instance of TProfile
	unsigned int binNumber = distributionIndices->GetBinNumber(0);
	string name = LinearisedDistribution->GetName();
	string title = LinearisedDistribution->GetTitle();
	TProfile * profilePlot = new TProfile( "temporaryProfile", "temporaryProfile", binNumber - 2, distributionIndices->GetBinLowEdgesForRoot( 0 ) );

	//Populate the profile
	vector< unsigned int > separateIndices;
	vector< double > centralValues;
	for ( unsigned int binIndex = 0; binIndex < distributionIndices->GetBinNumber(); binIndex++ )
	{
		//Work out the delinearised bin index and central value
		separateIndices = distributionIndices->GetNDimensionalIndex( binIndex );
		centralValues = distributionIndices->GetCentralValues( separateIndices );

		//Store the values
		profilePlot->Fill( centralValues[0], centralValues[1], LinearisedDistribution->GetBinContent( binIndex ) );
	}

	//Copy it into a histogram to avoid the stupid bug in Copy() with a profile cast as TH1D
	delete LinearisedDistribution;
	TH1F * returnHisto = new TH1F( name.c_str(), title.c_str(), binNumber - 2, distributionIndices->GetBinLowEdgesForRoot( 0 ) );
	for ( unsigned int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		double content = profilePlot->GetBinContent( binIndex );
		returnHisto->SetBinContent( binIndex, content );
	}

	//Return result
	delete profilePlot;
	return returnHisto;
}

vector< double > XvsYNormalisedPlotMaker::DelineariseErrors( vector< double > LinearisedErrors )
{
	//Find the target number of bins
	unsigned int binNumber = distributionIndices->GetBinNumber(0);

	//Make a vector of the de-linearised data
	vector< double > delinearisedErrors( binNumber, 0.0 );
	vector< unsigned int > separateIndices;
	for ( unsigned int binIndex = 0; binIndex < distributionIndices->GetBinNumber(); binIndex++ )
	{
		//Find the delinearised bin index
		separateIndices = distributionIndices->GetNDimensionalIndex(binIndex);

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
vector< double > XvsYNormalisedPlotMaker::CorrectedErrors()
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
TH2F * XvsYNormalisedPlotMaker::DAgostiniCovariance()
{
	if ( finalised )
	{
		return covarianceMatrix;
	}
	else
	{
		cerr << "Trying to retrieve D'Agostini covariance matrix from unfinalised XvsYNormalisedPlotMaker" << endl;
		exit(1);
	}
}

//Return the names of the variables involved
vector< string > XvsYNormalisedPlotMaker::VariableNames()
{
	vector< string > result( 1, xName );
	result.push_back( yName );
	return result;
}

//Return the type of correction the plot will perform
int XvsYNormalisedPlotMaker::CorrectionMode()
{
	return correctionType;
}

//Instantiate an object to correct the data
ICorrection * XvsYNormalisedPlotMaker::MakeCorrector( int CorrectionMode, IIndexCalculator * CorrectionIndices, string CorrectionName, unsigned int CorrectionIndex )
{
	if ( CorrectionMode == -1 )
	{
		return new Folding( CorrectionIndices, CorrectionName, CorrectionIndex );
	}
	else if ( CorrectionMode == 0 )
	{
		return new NoCorrection( CorrectionIndices, CorrectionName, CorrectionIndex );
	}
	else if ( CorrectionMode == 1 )
	{
		return new BinByBinUnfolding( CorrectionIndices, CorrectionName, CorrectionIndex );
	}
	else if ( CorrectionMode == 2 )
	{
		return new BayesianUnfolding( CorrectionIndices, CorrectionName, CorrectionIndex );
	}
	else
	{
		cerr << "Unrecognised correction mode (" << CorrectionMode << ") for " << CorrectionName << endl;
		exit(1);
	}
}

