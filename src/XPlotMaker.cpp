/**
  @class XPlotMaker

  Unfolds a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#include "XPlotMaker.h"
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
XPlotMaker::XPlotMaker()
{
}

//Constructor with the names to use for the variables
XPlotMaker::XPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum, int CorrectionMode, double ScaleFactor, bool Normalise )
{
	correctionType = CorrectionMode;
	xName = XVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	normalise = Normalise;
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

	//Make the x unfolder
	distributionIndices = new UniformIndices( binNumbers, minima, maxima );
	XUnfolder = MakeCorrector( correctionType );
}

//Constructor with the names to use for the variables
XPlotMaker::XPlotMaker( string XVariableName, string PriorName, vector< double > BinLowEdges, int CorrectionMode, double ScaleFactor, bool Normalise )
{
	correctionType = CorrectionMode;
	xName = XVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	normalise = Normalise;
	vector< vector< double > > binEdges;
	systematicRandom = 0;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Store the x range
	binEdges.push_back( BinLowEdges );

	//Make the x unfolder
	distributionIndices = new CustomIndices( binEdges );
	XUnfolder = MakeCorrector( correctionType );
}

//For use with Clone
XPlotMaker::XPlotMaker( string XVariableName, string PriorName, IIndexCalculator * DistributionIndices,
		unsigned int OriginalID, int CorrectionMode, double ScaleFactor, bool Normalise, vector< double > InputOffsets, vector< double > InputWidths )
{
	correctionType = CorrectionMode;
	xName = XVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	normalise = Normalise;

	//Set up for systematics
	systematicOffsets = InputOffsets;
	systematicWidths = InputWidths;
	if ( InputOffsets.size() == 0 )
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

	//Make the x unfolder
	distributionIndices = DistributionIndices;
	XUnfolder = MakeCorrector( correctionType );

	//Make the systematic unfolders too
	systematicUnfolders = vector< ICorrection* >( systematicOffsets.size(), XUnfolder->CloneShareSmearingMatrix() );
}

//Destructor
XPlotMaker::~XPlotMaker()
{
	delete distributionIndices;
	delete XUnfolder;
	for ( unsigned int experimentIndex = 0; experimentIndex < systematicUnfolders.size(); experimentIndex++ )
	{
		delete systematicUnfolders[ experimentIndex ];
	}
	if ( finalised )
	{
		delete correctedDistribution;
		delete uncorrectedDistribution;
		delete mcTruthDistribution;
		delete smearingMatrix;
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

//Set up a systematic error study
void XPlotMaker::AddSystematic( vector< double > SystematicOffset, vector< double > SystematicWidth, unsigned int NumberOfPseudoExperiments )
{
	//Check for correct number of observables
	if ( SystematicOffset.size() != 1 || SystematicWidth.size() != 1 )
	{
		cerr << "ERROR: There's only one observable in an XPlotmaker, so don't try sending multiple systematic offsets or widths at once" << endl;
		exit(1);
	}

	//Check for sensible pseudoexperiment number
	double experimentNumber = NumberOfPseudoExperiments;
	if ( SystematicWidth[0] == 0.0 && NumberOfPseudoExperiments > 1 )
	{
		cout << "INFO: There's no point doing many pseudoexperiments for a single systematic offset - using 1" << endl;
		experimentNumber = 1;
	}
	if ( SystematicWidth[0] != 0.0 && NumberOfPseudoExperiments < 1 )
	{
		cerr << "ERROR: You need to pick a number of pseudoexperiments to propagate this systematic width" << endl;
		exit(1);
	}

	//Make the random number generator
	if ( !systematicRandom && SystematicWidth[0] != 0.0 )
	{
		systematicRandom = new TRandom3(0);
	}

	//Make the extra experiments
	for ( unsigned int experimentIndex = 0; experimentIndex < experimentNumber; experimentIndex++ )
	{
		systematicUnfolders.push_back( XUnfolder->CloneShareSmearingMatrix() );
		systematicOffsets.push_back( SystematicOffset[0] );
		systematicWidths.push_back( SystematicWidth[0] );
	}
}

//Copy the object
XPlotMaker * XPlotMaker::Clone( string NewPriorName )
{
	return new XPlotMaker( xName, NewPriorName, distributionIndices->Clone(), thisPlotID, correctionType, scaleFactor, normalise, systematicOffsets, systematicWidths );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XPlotMaker::StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XPlotMaker" << endl;
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
		double truthWeight = TruthInput->EventWeight();
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x values
		truthValues.push_back( xTruthValue );
		reconstructedValues.push_back( xReconstructedValue );
		XUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );
	}
}
void XPlotMaker::StoreMiss( IFileInput * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XPlotMaker" << endl;
		exit(1);
	}
	else
	{
		vector< double > truthValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xTruthValue = TruthInput->GetValue( xName );
		double truthWeight = TruthInput->EventWeight();

		//Store the x value
		truthValues.push_back( xTruthValue );
		XUnfolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );
	}
}
void XPlotMaker::StoreFake( IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector< double > reconstructedValues;

		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( ReconstructedInput->Description() ) );

		//Retrieve the values from the Ntuple
		double xReconstructedValue = ReconstructedInput->GetValue( xName );
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x value
		reconstructedValues.push_back( xReconstructedValue );
		XUnfolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );
	}
}
void XPlotMaker::StoreData( IFileInput * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector< double > dataValues;
		double randomOne, randomTwo;

		//Retrieve the values from the Ntuple
		double xDataValue = DataInput->GetValue( xName );
		double dataWeight = DataInput->EventWeight();

		//Store the x value
		dataValues.push_back( xDataValue );
		XUnfolder->StoreDataValue( dataValues, dataWeight );

		//Do all the systematic error experiments
		for ( unsigned int experimentIndex = 0; experimentIndex < systematicOffsets.size(); experimentIndex++ )
		{
			dataValues[0] = xDataValue + systematicOffsets[ experimentIndex ];

			if ( systematicWidths[ experimentIndex ] != 0.0 )
			{
				systematicRandom->Rannor( randomOne, randomTwo );
				dataValues[0] += ( randomOne * systematicWidths[ experimentIndex ] );
			}

			systematicUnfolders[ experimentIndex ]->StoreDataValue( dataValues, dataWeight );
		}
	} 
}

//Do the unfolding
void XPlotMaker::Correct( unsigned int MostIterations, bool SkipUnfolding, unsigned int ErrorMode, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Unfold the distribution
		if ( !SkipUnfolding )
		{
			XUnfolder->Correct( MostIterations, ErrorMode, WithSmoothing );

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

		//Retrieve the results
		TH1F * XCorrected = XUnfolder->GetCorrectedHistogram( XFullName + "Corrected", XFullTitle + " Corrected Distribution", normalise );

		//And the systematics
		for ( unsigned int experimentIndex = 0; experimentIndex < systematicUnfolders.size(); experimentIndex++ )
		{
			stringstream systematicString;
			systematicString << thisPlotID << xName << priorName << "Systematic" << experimentIndex;

			systematicResults.push_back( systematicUnfolders[ experimentIndex ]->GetCorrectedHistogram( systematicString.str(), systematicString.str() ) );
		}

		//Retrieve some other bits for debug
		TH1F * XUncorrected = XUnfolder->GetUncorrectedHistogram( XFullName + "Uncorrected", XFullTitle + " Uncorrected Distribution", normalise );
		TH1F * XTruth = XUnfolder->GetTruthHistogram( XFullName + "Truth", XFullTitle + " Truth Distribution", normalise );
		TH2F * XSmearing = XUnfolder->GetSmearingHistogram( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );
		if ( ErrorMode > 1 && !SkipUnfolding )
		{
			covarianceMatrix = XUnfolder->DAgostiniCovariance( XFullName + "Covariance", XFullTitle + " Covariance Matrix" );
		}

		//Scale the histograms
		XCorrected->Scale( scaleFactor );
		XTruth->Scale( scaleFactor );
		XUncorrected->Scale( scaleFactor );

		//Get the error vectors
		vector< double > XErrors = XUnfolder->Variances();

		//Calculate errors
		TH1F * XCorrectedNotNormalised;
		if ( normalise )
		{
			XCorrectedNotNormalised = XUnfolder->GetCorrectedHistogram( XFullName + "CorrectedNotNormalised", XFullTitle + " Corrected Distribution Not Normalised", false );
			XCorrectedNotNormalised->Scale( scaleFactor );
		}
		for ( unsigned int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double combinedError = sqrt( XErrors[ binIndex ] ) * scaleFactor;

			//Scale the errors if the plots are normalised
			if ( normalise )
			{
				double normalisationFactor = XCorrected->GetBinContent( binIndex ) / XCorrectedNotNormalised->GetBinContent( binIndex );

				//Check for div0 errors
				if ( isinf( normalisationFactor ) )
				{
					normalisationFactor = 1.0;
				}
				if ( isnan( normalisationFactor ) )
				{
					normalisationFactor = 0.0;
				}

				combinedError *= normalisationFactor;
			}

			correctedDataErrors.push_back( combinedError );
		}

		//Format and save the corrected distribution
		correctedDistribution = XCorrected;
		correctedDistribution->SetXTitle( xName.c_str() );
		correctedDistribution->SetYTitle( "Events" );

		//Format and save the uncorrected distribution
		uncorrectedDistribution = XUncorrected;
		uncorrectedDistribution->SetXTitle( xName.c_str() );
		uncorrectedDistribution->SetYTitle( "Events" );

		//Format and save the truth distribution
		mcTruthDistribution = XTruth;
		mcTruthDistribution->SetXTitle( xName.c_str() );
		mcTruthDistribution->SetYTitle( "Events" );

		//Format and save the smearing matrix
		smearingMatrix = XSmearing;
		string smearingXTitle = xName + " Truth Bin";
		string smearingYTitle = xName + " Reconstructed Bin";
		smearingMatrix->SetXTitle( smearingXTitle.c_str() );
		smearingMatrix->SetYTitle( smearingYTitle.c_str() );

		if ( ErrorMode > 1 && !SkipUnfolding )
		{
			//Format the covariance matrix
			covarianceMatrix->SetXTitle( "Bin Number" );
			covarianceMatrix->SetYTitle( "Bin Number" );
		}

		//Bin-by-bin scaling of errors using the corrected data
		if ( correctionType != BAYESIAN_MODE || ErrorMode < 1 )
		{
			for ( unsigned int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
			{
				double errorScaleFactor = correctedDistribution->GetBinContent(binIndex) / uncorrectedDistribution->GetBinContent(binIndex);

				//Check for div0 errors
				if ( isinf( errorScaleFactor ) )
				{
					errorScaleFactor = 1.0;
				}
				if ( isnan( errorScaleFactor ) )
				{
					errorScaleFactor = 0.0;
				}

				correctedDataErrors[binIndex] *= errorScaleFactor;
			}
		}

		//Mark as done
		finalised = true;
	}
}

//Do a closure test
bool XPlotMaker::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
{
	return XUnfolder->ClosureTest( MostIterations, WithSmoothing );
}

//Make a cross-check with MC
unsigned int XPlotMaker::MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing )
{
	return XUnfolder->MonteCarloCrossCheck( InputPriorDistribution, InputSmearing, WithSmoothing );
}

//Return some plots
TH1F * XPlotMaker::CorrectedHistogram()
{
	if ( finalised )
	{
		return correctedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve corrected plot from unfinalised XPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XPlotMaker::UncorrectedHistogram()
{
	if ( finalised )
	{
		return uncorrectedDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve uncorrected plot from unfinalised XPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XPlotMaker::MCTruthHistogram()
{
	if ( finalised )
	{
		return mcTruthDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve MC truth plot from unfinalised XPlotMaker" << endl;
		exit(1);
	}
}

//Return a distribution for use in the cross-checks
Distribution * XPlotMaker::PriorDistributionForCrossCheck()
{
	return XUnfolder->GetTruthDistribution();
}
SmearingMatrix * XPlotMaker::SmearingMatrixForCrossCheck()
{
	return XUnfolder->GetSmearingMatrix();
}

TH2F * XPlotMaker::SmearingHistogram()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised XPlotMaker" << endl;
		exit(1);
	}
}

string XPlotMaker::Description( bool WithSpaces )
{
	return xName;
}
string XPlotMaker::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector< double > XPlotMaker::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedDataErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised XPlotMaker" << endl;
		exit(1);
	}
}
TH2F * XPlotMaker::DAgostiniCovariance()
{
	if ( finalised )
	{
		return covarianceMatrix;
	}
	else
	{
		cerr << "Trying to retrieve D'Agostini covariance matrix from unfinalised XPlotMaker" << endl;
		exit(1);
	}
}

//Return the names of the variables involved
vector< string > XPlotMaker::VariableNames()
{
	return vector< string >( 1, xName );
}

//Return the type of correction the plot will perform
int XPlotMaker::CorrectionMode()
{
	return correctionType;
}

//Return the results of the systematic experiments
vector< TH1F* > XPlotMaker::SystematicHistograms()
{
	return systematicResults;
}

//Instantiate an object to correct the data
ICorrection * XPlotMaker::MakeCorrector( int CorrectionMode )
{
	if ( CorrectionMode == -1 )
	{
		return new Folding( distributionIndices, xName + priorName, thisPlotID );
	}
	else if ( CorrectionMode == 0 )
	{
		return new NoCorrection( distributionIndices, xName + priorName, thisPlotID );
	}
	else if ( CorrectionMode == 1 )
	{
		return new BinByBinUnfolding( distributionIndices, xName + priorName, thisPlotID );
	}
	else if ( CorrectionMode == 2 )
	{
		return new BayesianUnfolding( distributionIndices, xName + priorName, thisPlotID );
	}
	else
	{
		cerr << "Unrecognised correction mode (" << CorrectionMode << ")" << endl;
		exit(1);
	}
}
