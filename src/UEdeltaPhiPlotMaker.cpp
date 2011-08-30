/**
  @class UEdeltaPhiPlotMaker

  Unfolds a 1D distribution which takes multiple values within a single event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 03-07-2011
 */

#include "UEdeltaPhiPlotMaker.h"
#include "BayesianUnfolding.h"
#include "UniformIndices.h"
#include "CustomIndices.h"
#include "Folding.h"
#include "NoCorrection.h"
#include "BinByBinUnfolding.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

const int BAYESIAN_MODE = 2;

//Default constructor - useless
UEdeltaPhiPlotMaker::UEdeltaPhiPlotMaker()
{
}

//Constructor with the names to use for the variables
UEdeltaPhiPlotMaker::UEdeltaPhiPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum,
		int CorrectionMode, vector< string > OtherVariableNames, double ScaleFactor, bool Normalise )
{
        cout << "WARNING: UEdeltaPhiPlotMaker is an experiment that I deem failed. Feel free to try it out, but I don't think it gives useful results" << endl;

	correctionType = CorrectionMode;
	xName = XVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	normalise = Normalise;
	vector< double > minima, maxima;
	vector< unsigned int > binNumbers;
	otherPairingNames = OtherVariableNames;
	otherPairingNames.push_back( xName );

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
	XUnfolder = MakeCorrector( CorrectionMode );
}

//Constructor with the names to use for the variables
UEdeltaPhiPlotMaker::UEdeltaPhiPlotMaker( string XVariableName, string PriorName, vector< double > BinLowEdges,
		int CorrectionMode, vector< string > OtherVariableNames, double ScaleFactor, bool Normalise )
{
        cout << "WARNING: UEdeltaPhiPlotMaker is an experiment that I deem failed. Feel free to try it out, but I don't think it gives useful results" << endl;

	correctionType = CorrectionMode;
	xName = XVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	normalise = Normalise;
	vector< vector< double > > binEdges;
	otherPairingNames = OtherVariableNames;
	otherPairingNames.push_back( xName );

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Store the x range
	binEdges.push_back( BinLowEdges );

	//Make the x unfolder
	distributionIndices = new CustomIndices( binEdges );
	XUnfolder = MakeCorrector( CorrectionMode );
}

//For use with Clone
UEdeltaPhiPlotMaker::UEdeltaPhiPlotMaker( vector< string > OtherVariableNames, string PriorName, IIndexCalculator * DistributionIndices,
		unsigned int OriginalID, int CorrectionMode, double ScaleFactor, bool Normalise )
{
        cout << "WARNING: UEdeltaPhiPlotMaker is an experiment that I deem failed. Feel free to try it out, but I don't think it gives useful results" << endl;

	correctionType = CorrectionMode;
	xName = OtherVariableNames[ OtherVariableNames.size() - 1 ];
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	normalise = Normalise;
	otherPairingNames = OtherVariableNames;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID + OriginalID;

	//Make the x unfolder
	distributionIndices = DistributionIndices;
	XUnfolder = MakeCorrector( CorrectionMode );
}

//Destructor
UEdeltaPhiPlotMaker::~UEdeltaPhiPlotMaker()
{
	delete distributionIndices;
	delete XUnfolder;
	otherPairingNames.clear();
	if ( finalised )
	{
		delete correctedDistribution;
		delete uncorrectedDistribution;
		delete mcTruthDistribution;
		delete smearingMatrix;
	}
}

//Copy the object
UEdeltaPhiPlotMaker * UEdeltaPhiPlotMaker::Clone( string NewPriorName )
{
	return new UEdeltaPhiPlotMaker( otherPairingNames, NewPriorName, distributionIndices->Clone(), thisPlotID, correctionType, scaleFactor, normalise );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void UEdeltaPhiPlotMaker::StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
	else
	{
		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Prepare the result vectors
		unpairedTruth.clear();
		unpairedReco.clear();
		truthRecoPairs.clear();
		vector< double > *truthValues, *reconstructedValues;

		//Loop over each supplied variable name to use in pairing, with the variable to plot as the last
		for ( unsigned int nameIndex = 0; nameIndex < otherPairingNames.size(); nameIndex++ )
		{
			//Retrieve the values from the Ntuple
			truthValues = TruthInput->GetVector( otherPairingNames[ nameIndex ] );
			reconstructedValues = ReconstructedInput->GetVector( otherPairingNames[ nameIndex ] );

			//Initialise the unpaired lists
			if ( nameIndex == 0 )
			{
				for ( unsigned int truthIndex = 0; truthIndex < truthValues->size(); truthIndex++ )
				{
					unpairedTruth.push_back( truthIndex );
				}
				for ( unsigned int reconstructedIndex = 0; reconstructedIndex < reconstructedValues->size(); reconstructedIndex++ )
				{
					unpairedReco.push_back( reconstructedIndex );
				}
			}

			//Do the pairing using this variable
			MakePairs( truthValues, reconstructedValues );
		}

		//Retrieve the weights from the Ntuple
		double truthWeight = TruthInput->EventWeight();
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Find the lead jet phi
		double truthLeadPhi = TruthInput->GetValue( "LeadJetPhi" );
		double reconstructedLeadPhi = ReconstructedInput->GetValue( "LeadJetPhi" );

		//Store the pairs
		for ( unsigned int pairIndex = 0; pairIndex < truthRecoPairs.size(); pairIndex++ )
		{
			double singleTruth = PlusMinusPi( ( *truthValues )[ truthRecoPairs[ pairIndex ].first ] - truthLeadPhi );
			double singleReco = PlusMinusPi( ( *reconstructedValues )[ truthRecoPairs[ pairIndex ].second ] - reconstructedLeadPhi );
			XUnfolder->StoreTruthRecoPair( vector< double >( 1, singleTruth ), vector< double >( 1, singleReco ), truthWeight, reconstructedWeight, useInPrior );
		}

		//Store the misses
		for ( unsigned int missIndex = 0; missIndex < unpairedTruth.size(); missIndex++ )
		{
			double singleTruth = PlusMinusPi( ( *truthValues )[ unpairedTruth[ missIndex ] ] - truthLeadPhi );
			XUnfolder->StoreUnreconstructedTruth( vector< double >( 1, singleTruth ), truthWeight, useInPrior );
		}

		//Store the fakes
		for ( unsigned int fakeIndex = 0; fakeIndex < unpairedReco.size(); fakeIndex++ )
		{
			double singleReco = PlusMinusPi( ( *reconstructedValues )[ unpairedReco[ fakeIndex ] ] - reconstructedLeadPhi );
			XUnfolder->StoreReconstructedFake( vector< double >( 1, singleReco ), reconstructedWeight, useInPrior );
		}
	}
}
void UEdeltaPhiPlotMaker::StoreMiss( IFileInput * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
	else
	{
		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Retrieve the values from the Ntuple
		vector< double > truthValues = *( TruthInput->GetVector( xName ) );
		double truthWeight = TruthInput->EventWeight();

		//Find the lead jet phi
		double truthLeadPhi = TruthInput->GetValue( "LeadJetPhi" );

		//Store the x values
		for ( unsigned int truthIndex = 0; truthIndex < truthValues.size(); truthIndex++ )
		{
			double singleTruth = PlusMinusPi( truthValues[ truthIndex ] - truthLeadPhi );
			XUnfolder->StoreUnreconstructedTruth( vector< double >( 1, singleTruth ), truthWeight, useInPrior );
		}
	}
}
void UEdeltaPhiPlotMaker::StoreFake( IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( ReconstructedInput->Description() ) );

		//Retrieve the values from the Ntuple
		vector< double > reconstructedValues = *( ReconstructedInput->GetVector( xName ) );
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Find the lead jet phi
		double reconstructedLeadPhi = ReconstructedInput->GetValue( "LeadJetPhi" );

		//Store the x value
		for ( unsigned int reconstructedIndex = 0; reconstructedIndex < reconstructedValues.size(); reconstructedIndex++ )
		{
			double singleReco = PlusMinusPi( reconstructedValues[ reconstructedIndex ] - reconstructedLeadPhi );
			XUnfolder->StoreReconstructedFake( vector< double >( 1, singleReco ), reconstructedWeight, useInPrior );
		}
	}
}
void UEdeltaPhiPlotMaker::StoreData( IFileInput * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		//Retrieve the values from the Ntuple
		vector< double > dataValues = *( DataInput->GetVector( xName ) );
		double dataWeight = DataInput->EventWeight();

		//Find the lead jet phi
		double dataLeadPhi = DataInput->GetValue( "LeadJetPhi" );

		//Store the x value
		for ( unsigned int dataIndex = 0; dataIndex < dataValues.size(); dataIndex++ )
		{
			double singleDatum = PlusMinusPi( dataValues[ dataIndex ] - dataLeadPhi );
			XUnfolder->StoreDataValue( vector< double >( 1, singleDatum ), dataWeight );
		}
	} 
}

//Do the unfolding
void UEdeltaPhiPlotMaker::Correct( unsigned int MostIterations, bool SkipUnfolding, unsigned int ErrorMode, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "UEdeltaPhiPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Unfold the distribution
		if ( !SkipUnfolding )
		{
			XUnfolder->Correct( MostIterations, ErrorMode, WithSmoothing );
		}

		//Make some plot titles
		stringstream uniqueIDString;
		uniqueIDString << thisPlotID;
		string XFullName = xName + priorName + uniqueIDString.str();
		string XFullTitle = xName + " using " + priorName;

		//Retrieve the results
		TH1F * XCorrected = XUnfolder->GetCorrectedHistogram( XFullName + "Corrected", XFullTitle + " Corrected Distribution", normalise );

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
bool UEdeltaPhiPlotMaker::ClosureTest( unsigned int MostIterations, bool WithSmoothing )
{
	return XUnfolder->ClosureTest( MostIterations, WithSmoothing );
}

//Make a cross-check with MC
unsigned int UEdeltaPhiPlotMaker::MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing )
{
	return XUnfolder->MonteCarloCrossCheck( InputPriorDistribution, InputSmearing, WithSmoothing );
}

//Return some plots
TH1F * UEdeltaPhiPlotMaker::CorrectedHistogram()
{
	if ( finalised )
	{
		return correctedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve corrected plot from unfinalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
}
TH1F * UEdeltaPhiPlotMaker::UncorrectedHistogram()
{
	if ( finalised )
	{
		return uncorrectedDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve uncorrected plot from unfinalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
}
TH1F * UEdeltaPhiPlotMaker::MCTruthHistogram()
{
	if ( finalised )
	{
		return mcTruthDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve MC truth plot from unfinalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
}

//Return a distribution for use in the cross-checks
Distribution * UEdeltaPhiPlotMaker::PriorDistributionForCrossCheck()
{
	return XUnfolder->GetTruthDistribution();
}
SmearingMatrix * UEdeltaPhiPlotMaker::SmearingMatrixForCrossCheck()
{
	return XUnfolder->GetSmearingMatrix();
}

TH2F * UEdeltaPhiPlotMaker::SmearingHistogram()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
}

string UEdeltaPhiPlotMaker::Description( bool WithSpaces )
{
	return xName;
}
string UEdeltaPhiPlotMaker::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector< double > UEdeltaPhiPlotMaker::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedDataErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
}
TH2F * UEdeltaPhiPlotMaker::DAgostiniCovariance()
{
	if ( finalised )
	{
		return covarianceMatrix;
	}
	else
	{
		cerr << "Trying to retrieve D'Agostini covariance matrix from unfinalised UEdeltaPhiPlotMaker" << endl;
		exit(1);
	}
}

//Return the names of the variables involved
vector< string > UEdeltaPhiPlotMaker::VariableNames()
{
	//return otherPairingNames;

	vector< string > relevantNames( otherPairingNames );
	relevantNames.push_back( "LeadJetPhi" );
	return relevantNames;
}

//Pair up values on the sub-event level using their proximity
void UEdeltaPhiPlotMaker::MakePairs( vector< double > * TruthValues, vector< double > * RecoValues )
{
	vector< vector< unsigned int > > guessPairings( RecoValues->size(), vector< unsigned int >() );

	for ( unsigned int truthIndex = 0; truthIndex < unpairedTruth.size(); truthIndex++ )
	{
		//Find the closest match 
		unsigned int closestRecoIndex;
		double smallestRecoDifference;
		for ( unsigned int recoIndex = 0; recoIndex < unpairedReco.size(); recoIndex++ )
		{
			double recoDifference = fabs( ( *TruthValues )[ unpairedTruth[ truthIndex ] ] - ( *RecoValues )[ unpairedReco[ recoIndex ] ] );

			if ( recoIndex == 0 || recoDifference < smallestRecoDifference )
			{
				closestRecoIndex = unpairedReco[ recoIndex ];
				smallestRecoDifference = recoDifference;
			}
		}

		//Store the match (if there was one to be made)
		if ( unpairedReco.size() > 0 )
		{
			guessPairings[ closestRecoIndex ].push_back( unpairedTruth[ truthIndex ] );
		}
	}

	//Look for truth-reco pairs
	for ( unsigned int recoIndex = 0; recoIndex < guessPairings.size(); recoIndex++ )
	{
		//Look for unconflicted pairs
		if ( guessPairings[ recoIndex ].size() > 0 )
		{
			//Find the truth that the reco most prefers
			unsigned int closestTruthIndex;
			double smallestTruthDifference;
			for ( unsigned int truthIndex = 0; truthIndex < guessPairings[ recoIndex ].size(); truthIndex++ )
			{
				double truthDifference = fabs( ( *TruthValues )[ guessPairings[ recoIndex ][ truthIndex ] ] - ( *RecoValues )[ recoIndex ] );

				if ( truthIndex == 0 || truthDifference < smallestTruthDifference )
				{
					closestTruthIndex = guessPairings[ recoIndex ][ truthIndex ];
					smallestTruthDifference = truthDifference;
				}
			}

			//Store the pair
			truthRecoPairs.push_back( make_pair( closestTruthIndex, recoIndex ) );

			//Remove the truth index from the unpaired truth list
			vector< unsigned int >::iterator searchIterator;
			for ( searchIterator = unpairedTruth.begin(); searchIterator != unpairedTruth.end(); searchIterator++ )
			{
				if ( *searchIterator == closestTruthIndex )
				{
					unpairedTruth.erase( searchIterator );
					break;
				}
			}

			//Remove the reco index from the unpaired reco list
			for ( searchIterator = unpairedReco.begin(); searchIterator != unpairedReco.end(); searchIterator++ )
			{
				if ( *searchIterator == recoIndex )
				{
					unpairedReco.erase( searchIterator );
					break;
				}
			}
		}
	}
}

double UEdeltaPhiPlotMaker::PlusMinusPi( double NotInRange )
{
	//Artificially limit the number of attempts to avoid floating point errors causing infinite loops
	for ( unsigned int attemptIndex = 0; attemptIndex < 4; attemptIndex++ )
	{
		if ( NotInRange > M_PI )
		{
			NotInRange -= 2.0 * M_PI;
		}
		else if ( NotInRange < -M_PI )
		{
			NotInRange += 2.0 * M_PI;
		}
		else
		{
			return NotInRange;
		}
	}

	return NotInRange;
}

//Return the type of correction the plot will perform
int UEdeltaPhiPlotMaker::CorrectionMode()
{
	return correctionType;
}

//Instantiate an object to correct the data
ICorrection * UEdeltaPhiPlotMaker::MakeCorrector( int CorrectionMode )
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
