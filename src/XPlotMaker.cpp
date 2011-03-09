/**
  @class XPlotMaker

  Unfolds a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#include "XPlotMaker.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "TFile.h"

using namespace std;

//Default constructor - useless
XPlotMaker::XPlotMaker()
{
}

//Constructor with the names to use for the variables
XPlotMaker::XPlotMaker( string XVariableName, string PriorName,
		int XBinNumber, double XMinimum, double XMaximum,
		double ScaleFactor, int UniqueID ) : xName( XVariableName ), priorName( PriorName ), finalised( false ), uniqueID(UniqueID), scaleFactor(ScaleFactor)
{
	vector<double> minima, maxima;
	vector<int> binNumbers;

	//Store the x range
	minima.push_back( XMinimum );
	maxima.push_back( XMaximum );
	binNumbers.push_back( XBinNumber );

	//Make the x unfolder
	XUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + priorName, uniqueID );

	//Set up the indices for the distributions
	DistributionIndices = new Indices( binNumbers, minima, maxima );
}

//Destructor
XPlotMaker::~XPlotMaker()
{
	delete DistributionIndices;
	delete XUnfolder;
	delete correctedDistribution;
	delete uncorrectedDistribution;
	delete mcTruthDistribution;
	delete smearingMatrix;
}

//Copy the object
IPlotMaker * XPlotMaker::Clone( string NewPriorName )
{
	return new XPlotMaker( xName, NewPriorName,
			DistributionIndices->GetBinNumber(0) - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0],
			scaleFactor, uniqueID );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XPlotMaker::StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XPlotMaker" << endl;
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
		double truthWeight = TruthInput->EventWeight();
		double reconstructedWeight = ReconstructedInput->EventWeight();

		//Store the x values
		truthValues.push_back( xTruthValue );
		reconstructedValues.push_back( xReconstructedValue );
		XUnfolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );
	}
}
void XPlotMaker::StoreMiss( InputNtuple * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XPlotMaker" << endl;
		exit(1);
	}
	else
	{
		vector<double> truthValues;

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
void XPlotMaker::StoreFake( InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector<double> reconstructedValues;

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
void XPlotMaker::StoreData( InputNtuple * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		vector<double> dataValues;

		//Retrieve the values from the Ntuple
		double xDataValue = DataInput->GetValue( xName );
		double dataWeight = DataInput->EventWeight();

		//Store the x value
		dataValues.push_back( xDataValue );
		XUnfolder->StoreDataValue( dataValues, dataWeight );
	} 
}

//Do the unfolding
void XPlotMaker::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Closure test
		XUnfolder->ClosureTest( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );

		//Unfold the distribution
		XUnfolder->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );

		//Make some plot titles
		char uniqueIDString[10];
		sprintf( uniqueIDString, "%d", uniqueID );
		string XFullName = xName + priorName + uniqueIDString;
		string XFullTitle = xName + " using " + priorName;

		//Retrieve the results
		TH1F * XCorrected = XUnfolder->GetUnfoldedHistogram( XFullName + "Corrected", XFullTitle + " Corrected Distribution" );

		//Retrieve some other bits for debug
		TH1F * XUncorrected = XUnfolder->GetUncorrectedDataHistogram( XFullName + "Uncorrected", XFullTitle + " Uncorrected Distribution" );
		TH1F * XTruth = XUnfolder->GetTruthHistogram( XFullName + "Truth", XFullTitle + " Truth Distribution" );
		TH2F * XSmearing = XUnfolder->GetSmearingMatrix( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );

		//Scale the histograms
		XCorrected->Scale( scaleFactor );
		XTruth->Scale( scaleFactor );
		XUncorrected->Scale( scaleFactor );

		//Get the error vectors
		vector<double> XErrors = XUnfolder->SumOfDataWeightSquares();

		//Calculate errors
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double combinedError = sqrt( XErrors[ binIndex ] ) * scaleFactor;
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
		string smearingXTitle = xName + " Truth";
		string smearingYTitle = xName + " Reconstructed";
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
int XPlotMaker::MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	return XUnfolder->MonteCarloCrossCheck( ReferenceDistribution, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
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
Distribution * XPlotMaker::MonteCarloTruthForCrossCheck()
{
	return XUnfolder->GetTruthDistribution();
}

TH2F * XPlotMaker::SmearingMatrix()
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
vector<double> XPlotMaker::CorrectedErrors()
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

