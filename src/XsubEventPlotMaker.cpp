/**
  @class XsubEventPlotMaker

  Unfolds a 1D distribution which takes multiple values within a single event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 03-07-2011
 */

#include "XsubEventPlotMaker.h"
#include "UniformIndices.h"
#include "CustomIndices.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

//Default constructor - useless
XsubEventPlotMaker::XsubEventPlotMaker()
{
}

//Constructor with the names to use for the variables
XsubEventPlotMaker::XsubEventPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum, vector< string > OtherVariableNames, double ScaleFactor, bool Normalise )
{
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
	XUnfolder = new IterativeUnfolding( distributionIndices, xName + priorName, thisPlotID );
}

//Constructor with the names to use for the variables
XsubEventPlotMaker::XsubEventPlotMaker( string XVariableName, string PriorName, vector< double > BinLowEdges, vector< string > OtherVariableNames, double ScaleFactor, bool Normalise )
{
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
	XUnfolder = new IterativeUnfolding( distributionIndices, xName + priorName, thisPlotID );
}

//For use with Clone
XsubEventPlotMaker::XsubEventPlotMaker( vector< string > OtherVariableNames, string PriorName, IIndexCalculator * DistributionIndices,
		unsigned int OriginalID, double ScaleFactor, bool Normalise )
{
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
	XUnfolder = new IterativeUnfolding( distributionIndices, xName + priorName, thisPlotID );
}

//Destructor
XsubEventPlotMaker::~XsubEventPlotMaker()
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
IUnfolder * XsubEventPlotMaker::Clone( string NewPriorName )
{
	return new XsubEventPlotMaker( otherPairingNames, NewPriorName, distributionIndices->Clone(), thisPlotID, scaleFactor, normalise );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XsubEventPlotMaker::StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XsubEventPlotMaker" << endl;
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

		//Store the pairs
		for ( unsigned int pairIndex = 0; pairIndex < truthRecoPairs.size(); pairIndex++ )
		{
			double singleTruth = ( *truthValues )[ truthRecoPairs[ pairIndex ].first ];
			double singleReco = ( *reconstructedValues )[ truthRecoPairs[ pairIndex ].second ];
			XUnfolder->StoreTruthRecoPair( vector< double >( 1, singleTruth ), vector< double >( 1, singleReco ), truthWeight, reconstructedWeight, useInPrior );
		}

		//Store the misses
		for ( unsigned int missIndex = 0; missIndex < unpairedTruth.size(); missIndex++ )
		{
			double singleTruth = ( *truthValues )[ unpairedTruth[ missIndex ] ];
			XUnfolder->StoreUnreconstructedTruth( vector< double >( 1, singleTruth ), truthWeight, useInPrior );
		}

		//Store the fakes
		for ( unsigned int fakeIndex = 0; fakeIndex < unpairedReco.size(); fakeIndex++ )
		{
			double singleReco = ( *reconstructedValues )[ unpairedReco[ fakeIndex ] ];
			XUnfolder->StoreReconstructedFake( vector< double >( 1, singleReco ), reconstructedWeight, useInPrior );
		}
	}
}
void XsubEventPlotMaker::StoreMiss( IFileInput * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XsubEventPlotMaker" << endl;
		exit(1);
	}
	else
	{
		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( TruthInput->Description() ) );

		//Retrieve the values from the Ntuple
		vector< double > truthValues = *( TruthInput->GetVector( xName ) );
		double truthWeight = TruthInput->EventWeight();

		//Store the x values
		for ( unsigned int truthIndex = 0; truthIndex < truthValues.size(); truthIndex++ )
		{
			XUnfolder->StoreUnreconstructedTruth( vector< double >( 1, truthValues[ truthIndex ] ), truthWeight, useInPrior );
		}
	}
}
void XsubEventPlotMaker::StoreFake( IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XsubEventPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		//Find out if this is the correct prior
		bool useInPrior = ( priorName == *( ReconstructedInput->Description() ) );

		//Retrieve the values from the Ntuple
		vector< double > reconstructedValues = *( ReconstructedInput->GetVector( xName ) );
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x value
		for ( unsigned int reconstructedIndex = 0; reconstructedIndex < reconstructedValues.size(); reconstructedIndex++ )
		{
			XUnfolder->StoreReconstructedFake( vector< double >( 1, reconstructedValues[ reconstructedIndex ] ), reconstructedWeight, useInPrior );
		}
	}
}
void XsubEventPlotMaker::StoreData( IFileInput * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XsubEventPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		//Retrieve the values from the Ntuple
		vector< double > dataValues = *( DataInput->GetVector( xName ) );
		double dataWeight = DataInput->EventWeight();

		//Store the x value
		for ( unsigned int dataIndex = 0; dataIndex < dataValues.size(); dataIndex++ )
		{
			XUnfolder->StoreDataValue( vector< double >( 1, dataValues[ dataIndex ] ), dataWeight );
		}
	} 
}

//Do the unfolding
void XsubEventPlotMaker::Unfold( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool SkipUnfolding, unsigned int ErrorMode, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XsubEventPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Unfold the distribution
		if ( !SkipUnfolding )
		{
			XUnfolder->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, ErrorMode, WithSmoothing );
		}

		//Make some plot titles
		stringstream uniqueIDString;
		uniqueIDString << thisPlotID;
		string XFullName = xName + priorName + uniqueIDString.str();
		string XFullTitle = xName + " using " + priorName;

		//Retrieve the results
		TH1F * XCorrected = XUnfolder->GetUnfoldedHistogram( XFullName + "Corrected", XFullTitle + " Corrected Distribution", normalise );

		//Retrieve some other bits for debug
		TH1F * XUncorrected = XUnfolder->GetUncorrectedDataHistogram( XFullName + "Uncorrected", XFullTitle + " Uncorrected Distribution", normalise );
		TH1F * XTruth = XUnfolder->GetTruthHistogram( XFullName + "Truth", XFullTitle + " Truth Distribution", normalise );
		TH2F * XSmearing = XUnfolder->GetSmearingMatrix( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );
		if ( ErrorMode > 1 && !SkipUnfolding )
		{
			covarianceMatrix = XUnfolder->DAgostiniCovariance( XFullName + "Covariance", XFullTitle + " Covariance Matrix" );
		}

		//Scale the histograms
		XCorrected->Scale( scaleFactor );
		XTruth->Scale( scaleFactor );
		XUncorrected->Scale( scaleFactor );

		//Get the error vectors
		vector< double > XErrors = XUnfolder->SumOfDataWeightSquares();
		vector< double > XVariance;
		if ( ErrorMode > 0 && !SkipUnfolding )
		{
			XVariance = XUnfolder->DAgostiniVariance();
		}

		//Calculate errors
		TH1F * XCorrectedNotNormalised;
		if ( normalise )
		{
			XCorrectedNotNormalised = XUnfolder->GetUnfoldedHistogram( XFullName + "CorrectedNotNormalised", XFullTitle + " Corrected Distribution Not Normalised", false );
			XCorrectedNotNormalised->Scale( scaleFactor );
		}
		for ( unsigned int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double combinedError = sqrt( XErrors[ binIndex ] ) * scaleFactor;

			double dagostiniError;
			if ( ErrorMode > 0 && !SkipUnfolding )
			{
				dagostiniError = sqrt( XVariance[ binIndex ] ) * scaleFactor;
			}

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

				if ( ErrorMode > 0 && !SkipUnfolding )
				{
					dagostiniError *= normalisationFactor;
				}
			}

			correctedDataErrors.push_back( combinedError );

			if ( ErrorMode > 0 && !SkipUnfolding )
			{
				dagostiniErrors.push_back( dagostiniError );
			}
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

		//Mark as done
		finalised = true;
	}
}

//Do a closure test
bool XsubEventPlotMaker::ClosureTest( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	return XUnfolder->ClosureTest( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
}

//Make a cross-check with MC
unsigned int XsubEventPlotMaker::MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	return XUnfolder->MonteCarloCrossCheck( ReferenceDistribution, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
}

//Return some plots
TH1F * XsubEventPlotMaker::CorrectedHistogram()
{
	if ( finalised )
	{
		return correctedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve corrected plot from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XsubEventPlotMaker::UncorrectedHistogram()
{
	if ( finalised )
	{
		return uncorrectedDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve uncorrected plot from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}
TH1F * XsubEventPlotMaker::MCTruthHistogram()
{
	if ( finalised )
	{
		return mcTruthDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve MC truth plot from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}

//Return a distribution for use in the cross-checks
Distribution * XsubEventPlotMaker::MonteCarloTruthForCrossCheck()
{
	return XUnfolder->GetTruthDistribution();
}

TH2F * XsubEventPlotMaker::SmearingMatrix()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}

string XsubEventPlotMaker::Description( bool WithSpaces )
{
	return xName;
}
string XsubEventPlotMaker::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector< double > XsubEventPlotMaker::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedDataErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}
vector< double > XsubEventPlotMaker::DAgostiniErrors()
{
	if ( finalised )
	{
		return dagostiniErrors;
	}
	else
	{
		cerr << "Trying to retrieve D'Agostini errors from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}
TH2F * XsubEventPlotMaker::DAgostiniCovariance()
{
	if ( finalised )
	{
		return covarianceMatrix;
	}
	else
	{
		cerr << "Trying to retrieve D'Agostini covariance matrix from unfinalised XsubEventPlotMaker" << endl;
		exit(1);
	}
}

//Return the names of the variables involved
vector< string > XsubEventPlotMaker::VariableNames()
{
	return otherPairingNames;
}

//Pair up values on the sub-event level using their proximity
void XsubEventPlotMaker::MakePairs( vector< double > * TruthValues, vector< double > * RecoValues )
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

		//Store the match
		if ( closestRecoIndex >= guessPairings.size() )
		{
			cerr << closestRecoIndex << " >= " << guessPairings.size() << endl;
		}
		if ( unpairedReco.size() > 0 )
		{
			guessPairings[ closestRecoIndex ].push_back( unpairedTruth[ truthIndex ] );
		}
	}

	//Look for truth-reco pairs
	for ( unsigned int recoIndex = 0; recoIndex < guessPairings.size(); recoIndex++ )
	{
		//Look for unconflicted pairs
		if ( guessPairings[ recoIndex ].size() == 1 )
		{
			//Store the pair
			truthRecoPairs.push_back( make_pair( guessPairings[ recoIndex ][ 0 ], recoIndex ) );

			//Remove the truth index from the unpaired truth list
			vector< unsigned int >::iterator searchIterator;
			for ( searchIterator = unpairedTruth.begin(); searchIterator != unpairedTruth.end(); searchIterator++ )
			{
				if ( *searchIterator == guessPairings[ recoIndex ][ 0 ] )
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
