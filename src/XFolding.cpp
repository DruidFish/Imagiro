/**
  @class XFolding

  Applies the smearing matrix to a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#include "XFolding.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

//Default constructor - useless
XFolding::XFolding()
{
}

//Constructor with the names to use for the variables
XFolding::XFolding( string XVariableName, string PriorName,
		int XBinNumber, double XMinimum, double XMaximum,
		double ScaleFactor ) : xName( XVariableName ), priorName( PriorName ), finalised( false ), scaleFactor(ScaleFactor)
{
	vector<double> minima, maxima;
	vector<int> binNumbers;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Store the x range
	minima.push_back( XMinimum );
	maxima.push_back( XMaximum );
	binNumbers.push_back( XBinNumber );

	//Make the x unfolder
	XFolder = new Folding( binNumbers, minima, maxima, xName + priorName, thisPlotID );

	//Set up the indices for the distributions
	DistributionIndices = new Indices( binNumbers, minima, maxima );
}

//Destructor
XFolding::~XFolding()
{
	delete DistributionIndices;
	delete XFolder;
	delete foldedDistribution;
	delete inputDistribution;
	delete reconstructedDistribution;
	delete smearingMatrix;
}

//Copy the object
IPlotMaker * XFolding::Clone( string NewPriorName )
{
	return new XFolding( xName, NewPriorName,
			DistributionIndices->GetBinNumber(0) - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0],
			scaleFactor );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XFolding::StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XFolding" << endl;
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
		XFolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );
	}
}
void XFolding::StoreMiss( InputNtuple * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XFolding" << endl;
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
		XFolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );
	}
}
void XFolding::StoreFake( InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XFolding" << endl;
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
		XFolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );
	}
}
void XFolding::StoreData( InputNtuple * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XFolding" << endl;
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
		XFolder->StoreValueToFold( dataValues, dataWeight );
	} 
}

//Do the unfolding
void XFolding::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "XFolding is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Closure tests
		XFolder->ClosureTest();

		//Unfold the distributions
		XFolder->Fold();

		//Make some plot titles
		stringstream uniqueIDString;
		uniqueIDString << thisPlotID;
		string XFullName = xName + priorName + uniqueIDString.str();
		string XFullTitle = xName + " using " + priorName;

		//Retrieve the results
		TH1F * XSmeared = XFolder->GetFoldedHistogram( XFullName + "Smeared", XFullTitle + " Smeared Distribution" );

		//Retrieve some other bits for debug
		TH1F * XNotSmeared = XFolder->GetInputHistogram( XFullName + "Input", XFullTitle + " Input Distribution" );
		TH1F * XTruth = XFolder->GetReconstructedHistogram( XFullName + "Reconstructed", XFullTitle + " Reconstructed Distribution" );
		TH2F * XSmearing = XFolder->GetSmearingMatrix( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );

		//Scale the histograms
		XSmeared->Scale( scaleFactor );
		XTruth->Scale( scaleFactor );
		XNotSmeared->Scale( scaleFactor );

		//Get the error vector
		vector<double> XErrors = XFolder->SumOfInputWeightSquares();

		//Combine errors
		for ( int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double combinedError = sqrt( XErrors[ binIndex ] ) * scaleFactor;
			correctedInputErrors.push_back( combinedError );
		}

		//Format and save the corrected distribution
		foldedDistribution = XSmeared;
		foldedDistribution->SetXTitle( xName.c_str() );
		foldedDistribution->SetYTitle( "Events" );

		//Format and save the uncorrected distribution
		inputDistribution = XNotSmeared;
		inputDistribution->SetXTitle( xName.c_str() );
		inputDistribution->SetYTitle( "Events" );

		//Format and save the truth distribution
		reconstructedDistribution = XTruth;
		reconstructedDistribution->SetXTitle( xName.c_str() );
		reconstructedDistribution->SetXTitle( "Events" );

		//Format and save the smearing matrix
		smearingMatrix = XSmearing;
		string smearingXTitle = xName + " Truth";
		string smearingYTitle = xName + " Reconstructed";
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
int XFolding::MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	//Meaningless values
	ChiSquaredThreshold = 10.0;
	KolmogorovThreshold = 0.1;
	return 10;
}

//Return some plots
TH1F * XFolding::CorrectedHistogram()
{
	if ( finalised )
	{
		return foldedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve smeared plot from unfinalised XFolding" << endl;
		exit(1);
	}
}
TH1F * XFolding::UncorrectedHistogram()
{
	if ( finalised )
	{
		return inputDistribution;
	}       
	else
	{
		cerr << "Trying to retrieve plot without smearing from unfinalised XFolding" << endl;
		exit(1);
	}
}

TH1F * XFolding::MCTruthHistogram()
{
	if ( finalised )
	{
		return reconstructedDistribution;
	}
	else
	{
		cerr << "Trying to retrieve MC truth plot from unfinalised XFolding" << endl;
		exit(1);
	}
}

//Return a distribution for use in the cross-checks
Distribution * XFolding::MonteCarloTruthForCrossCheck() 
{
	//This doesn't really make sense to use
	return XFolder->GetTruthDistribution();
}


TH2F * XFolding::SmearingMatrix()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised XFolding" << endl;
		exit(1);
	}
}

string XFolding::Description( bool WithSpaces )
{
	return xName;
}
string XFolding::PriorName()
{
	return priorName;
}

//Error info for corrected distribution
vector<double> XFolding::CorrectedErrors()
{
	if ( finalised )
	{
		return correctedInputErrors;
	}
	else
	{
		cerr << "Trying to retrieve corrected data errors from unfinalised XFolding" << endl;
		exit(1);
	}
}

