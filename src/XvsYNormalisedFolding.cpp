/**
  @class XvsYNormalisedFolding

  Folds a 2D distribution, and divides it by the folded 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#include "XvsYNormalisedFolding.h"
#include "UniformIndices.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

//Default constructor - useless
XvsYNormalisedFolding::XvsYNormalisedFolding()
{
}

//Constructor with the names to use for the variables
XvsYNormalisedFolding::XvsYNormalisedFolding( string XVariableName, string YVariableName, string PriorName,
		unsigned int XBinNumber, double XMinimum, double XMaximum,
		unsigned int YBinNumber, double YMinimum, double YMaximum, double ScaleFactor )
{
	xName = XVariableName;
	yName = YVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;
	vector< double > minima, maxima;
	vector< unsigned int > binNumbers;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Store the x range
	minima.push_back( XMinimum );
	maxima.push_back( XMaximum );
	binNumbers.push_back( XBinNumber );

	//Make the x unfolder
	xIndices = new UniformIndices( binNumbers, minima, maxima );
	XFolder = new Folding( xIndices, xName + priorName, thisPlotID );

	//Store the y range
	minima.push_back( YMinimum );
	maxima.push_back( YMaximum );
	binNumbers.push_back( YBinNumber );

	//Set up the indices for the distributions
	distributionIndices = new UniformIndices( binNumbers, minima, maxima );

	//Make the x vs y unfolder
	XvsYFolder = new Folding( distributionIndices, xName + "vs" + yName + priorName, thisPlotID );

	//Set up the cross-check for data loss in delinearisation
	stringstream idString;
	idString << thisPlotID;
	string xvsyReconstructionName = xName + yName + priorName + "ReconstructionCheck" + idString.str();
	string xReconstructionName = xName + priorName + "ReconstructionCheck" + idString.str();
	xvsyReconstructionCheck = new TH1F( xvsyReconstructionName.c_str(), xvsyReconstructionName.c_str(), XBinNumber, XMinimum, XMaximum );
	xReconstructionCheck = new TH1F( xReconstructionName.c_str(), xReconstructionName.c_str(), XBinNumber, XMinimum, XMaximum );

	//Make a summary of the y input values
	yValueSummary = new StatisticsSummary();
}

XvsYNormalisedFolding::XvsYNormalisedFolding( string XVariableName, string YVariableName, string PriorName, IIndexCalculator * XIndices, IIndexCalculator * DistributionIndices, double ScaleFactor )
{
	xName = XVariableName;
	yName = YVariableName;
	priorName = PriorName;
	finalised = false;
	scaleFactor = ScaleFactor;

	//Set up a variable to keep track of the number of plots - used to prevent Root from complaining about making objects with the same names
	static unsigned int uniqueID = 0;
	uniqueID++;
	thisPlotID = uniqueID;

	//Set up the indices for the distributions
	xIndices = XIndices;
	distributionIndices = DistributionIndices;

	//Make the x unfolder
	XFolder = new Folding( xIndices, xName + priorName, thisPlotID );

	//Make the x vs y unfolder
	XvsYFolder = new Folding( distributionIndices, xName + "vs" + yName + priorName, thisPlotID );

	//Set up the cross-check for data loss in delinearisation
	stringstream idString;
	idString << thisPlotID;
	string xvsyReconstructionName = xName + yName + priorName + "ReconstructionCheck" + idString.str();
	string xReconstructionName = xName + priorName + "ReconstructionCheck" + idString.str();
	xvsyReconstructionCheck = new TH1F( xvsyReconstructionName.c_str(), xvsyReconstructionName.c_str(), distributionIndices->GetBinNumber(0) - 2, distributionIndices->GetBinLowEdgesForRoot(0) );
	xReconstructionCheck = new TH1F( xReconstructionName.c_str(), xReconstructionName.c_str(), xIndices->GetBinNumber(0) - 2, xIndices->GetBinLowEdgesForRoot(0) );

	//Make a summary of the y input values
	yValueSummary = new StatisticsSummary();
}

//Destructor
XvsYNormalisedFolding::~XvsYNormalisedFolding()
{
	delete xIndices;
	delete distributionIndices;
	delete XFolder;
	delete XvsYFolder;
	delete foldedDistribution;
	delete inputDistribution;
	delete reconstructedDistribution;
	delete smearingMatrix;
}

//Copy the object
IFolder * XvsYNormalisedFolding::Clone( string NewPriorName )
{
	return new XvsYNormalisedFolding( xName, yName, NewPriorName, xIndices->Clone(), distributionIndices->Clone(), scaleFactor );
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void XvsYNormalisedFolding::StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised XvsYNormalisedFolding" << endl;
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
		XFolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		//Store the y values
		truthValues.push_back( yTruthValue );
		reconstructedValues.push_back( yReconstructedValue );
		XvsYFolder->StoreTruthRecoPair( truthValues, reconstructedValues, truthWeight, reconstructedWeight, useInPrior );

		if ( useInPrior )
		{
			//Store values for checking the delinearisation
			xReconstructionCheck->Fill( xReconstructedValue, reconstructedWeight );
			xvsyReconstructionCheck->Fill( xReconstructedValue, yReconstructedValue * reconstructedWeight );
		}
	}
}
void XvsYNormalisedFolding::StoreMiss( IFileInput * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised XvsYNormalisedFolding" << endl;
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
		XFolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );

		//Store the y value
		truthValues.push_back( yTruthValue );
		XvsYFolder->StoreUnreconstructedTruth( truthValues, truthWeight, useInPrior );
	}
}
void XvsYNormalisedFolding::StoreFake( IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised XvsYNormalisedFolding" << endl;
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
		XFolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );

		//Store the y value
		reconstructedValues.push_back( yReconstructedValue );
		XvsYFolder->StoreReconstructedFake( reconstructedValues, reconstructedWeight, useInPrior );

		if ( useInPrior )
		{
			//Store values for checking the delinearisation
			xReconstructionCheck->Fill( xReconstructedValue, reconstructedWeight );
			xvsyReconstructionCheck->Fill( xReconstructedValue, yReconstructedValue * reconstructedWeight );
		}
	}
}
void XvsYNormalisedFolding::StoreData( IFileInput * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised XvsYNormalisedFolding" << endl;
		exit(1);
	}       
	else
	{
		vector< double > dataValues;

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
		yValueSummary->StoreEvent( yDataValue, dataWeight );

		//Store values for performing the delinearisation
		distributionIndices->StoreDataValue( dataValues, dataWeight );
	} 
}

//Do the unfolding
void XvsYNormalisedFolding::Fold()
{
	if ( finalised )
	{
		cerr << "XvsYNormalisedFolding is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Unfold the distributions
		XFolder->Fold();
		XvsYFolder->Fold();

		//Make some plot titles
		stringstream uniqueIDString;
		uniqueIDString << thisPlotID;
		string XFullName = xName + priorName + uniqueIDString.str();
		string XFullTitle = xName + " using " + priorName;
		string XvsYName = xName + "vs" + yName + priorName + uniqueIDString.str();
		string XvsYTitle = xName + " vs " + yName + " using " + priorName;

		//Retrieve the results
		TH1F * XSmeared = XFolder->GetFoldedHistogram( XFullName + "Smeared", XFullTitle + " Smeared Distribution" );
		TH1F * XvsYSmeared = XvsYFolder->GetFoldedHistogram( XvsYName + "Smeared", XvsYTitle + " Smeared Distribution" );

		//Retrieve some other bits for debug
		TH1F * XNotSmeared = XFolder->GetInputHistogram( XFullName + "Input", XFullTitle + " Input Distribution" );
		TH1F * XvsYNotSmeared = XvsYFolder->GetInputHistogram( XvsYName + "Input", XvsYTitle + " Input Distribution" );

		TH1F * XTruth = XFolder->GetReconstructedHistogram( XFullName + "Reconstructed", XFullTitle + " Reconstructed Distribution" );
		TH1F * XvsYTruth = XvsYFolder->GetReconstructedHistogram( XvsYName + "Reconstructed", XvsYTitle + " Reconstructed Distribution" );

		TH2F * XvsYSmearing = XvsYFolder->GetSmearingMatrix( XvsYName + "Smearing", XvsYTitle + " Smearing Matrix" );

		//De-linearise the x vs y distributions
		TH1F * DelinearisedXvsYSmeared = Delinearise(XvsYSmeared);
		TH1F * DelinearisedXvsYNotSmeared = Delinearise(XvsYNotSmeared);
		TH1F * DelinearisedXvsYTruth = Delinearise(XvsYTruth);

		//Get the error vectors
		vector< double > XErrors = XFolder->SumOfInputWeightSquares();
		vector< double > XvsYErrors = XvsYFolder->SumOfInputWeightSquares();

		//Delinearise the x vs y errors
		vector< double > delinearisedXvsYErrors = DelineariseErrors( XvsYErrors );

		//Combine errors
		for ( unsigned int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
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

		//Check for data loss in the delinearisation
		xvsyReconstructionCheck->Divide( xvsyReconstructionCheck, xReconstructionCheck, scaleFactor, 1.0 );
		bool dataLost = false;
		double averagePercentError = 0.0;
		for ( int binIndex = 0; binIndex < xReconstructionCheck->GetNbinsX() + 2; binIndex++ )
		{
			//Compare the delinearised value with one calculated without going through delinearisation
			double correctValue = xvsyReconstructionCheck->GetBinContent( binIndex );
			double delinearisedValue = DelinearisedXvsYTruth->GetBinContent( binIndex );
			double errorFraction = fabs( delinearisedValue - correctValue ) / correctValue;
			double percentError = errorFraction * 100.0;

			//Check for stupid values
			if ( isnan( percentError ) )
			{
				percentError = 0.0;
			}

			//Add this error to the overall error calculation
			correctedInputErrors[ binIndex ] *= ( 1.0 + errorFraction );

			//Increment average error
			averagePercentError += percentError;

			//If there's a > 1%  discrepancy between the delinearised truth value and the correct value, warn that data is being lost in binning
			if ( fabs( percentError ) > 1.0 )
			{
				cerr << "Bin " << binIndex << " has value " << delinearisedValue << " which has a " << percentError << "\% error vs the reference value (" << correctValue << ")" << endl;
				dataLost = true;
			}
		}
		averagePercentError /= (double)( xReconstructionCheck->GetNbinsX() + 2 );
		cout << "Average bin error from delinearisation: " << averagePercentError << "\%" << endl;
		if ( dataLost || averagePercentError > 0.5 )
		{
			cerr << "Delinearisation has caused significant changes in distribution bins compared to their reference values" << endl;
			cerr << "Improve the binning of " << yName << " in the " << xName << " vs " << yName << " plot" << endl;
			cerr << "Try using a finer binning, and avoid large numbers of events in under or overflow bins" << endl;
			cerr << "Suggested bin width: " << yValueSummary->OptimumBinWidth() << endl;
			exit(1);
		}

		//Free some memory
		delete XSmeared;
		delete XNotSmeared;
		delete XTruth;
		delete yValueSummary;

		//Get the y range to plot
		//double yMinimum = distributionIndices->GetMinima()[1] * scaleFactor;
		//double yMaximum = distributionIndices->GetMaxima()[1] * scaleFactor;

		//Format and save the corrected distribution
		foldedDistribution = DelinearisedXvsYSmeared;
		foldedDistribution->SetXTitle( xName.c_str() );
		foldedDistribution->SetYTitle( yName.c_str() );
		//foldedDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the uncorrected distribution
		inputDistribution = DelinearisedXvsYNotSmeared;
		inputDistribution->SetXTitle( xName.c_str() );
		inputDistribution->SetYTitle( yName.c_str() );
		//inputDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the truth distribution
		reconstructedDistribution = DelinearisedXvsYTruth;
		reconstructedDistribution->SetXTitle( xName.c_str() );
		reconstructedDistribution->SetYTitle( yName.c_str() );
		//reconstructedDistribution->GetYaxis()->SetRangeUser( yMinimum, yMaximum );

		//Format and save the smearing matrix
		smearingMatrix = XvsYSmearing;
		string smearingXTitle = xName + " vs " + yName + " Truth";
		string smearingYTitle = xName + " vs " + yName + " Reconstructed";
		smearingMatrix->SetXTitle( smearingXTitle.c_str() );
		smearingMatrix->SetYTitle( smearingYTitle.c_str() );

		//Bin-by-bin scaling of errors using the corrected data
		for ( unsigned int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
		{
			double errorScaleFactor = foldedDistribution->GetBinContent(binIndex) / inputDistribution->GetBinContent(binIndex);

			//Check for div0 errors
			if ( isinf( errorScaleFactor ) )
			{
				errorScaleFactor = 1.0;
			}
			if ( isnan( errorScaleFactor ) )
			{
				errorScaleFactor = 0.0;
			}

			correctedInputErrors[binIndex] *= errorScaleFactor;
		}

		//Mark as done
		finalised = true;
	}
}

//Do a closure test
bool XvsYNormalisedFolding::ClosureTest()
{
	bool result = XvsYFolder->ClosureTest();
	result &= XFolder->ClosureTest();
	return result;
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
	unsigned int binNumber = distributionIndices->GetBinNumber(0);

	//Make a vector of the de-linearised data
	vector< double > delinearisedDistribution( binNumber, 0.0 );
	vector< unsigned int > separateIndices;
	vector< double > centralValues, dataCentralValues;
	for ( unsigned int binIndex = 0; binIndex < distributionIndices->GetBinNumber(); binIndex++ )
	{
		//Work out the delinearised bin index and central value
		separateIndices = distributionIndices->GetNDimensionalIndex(binIndex);
		centralValues = distributionIndices->GetCentralValues(separateIndices);

		//Increment the bin in the delinearised distribution
		double thisBinValue = LinearisedDistribution->GetBinContent(binIndex) * centralValues[1];
		delinearisedDistribution[ separateIndices[0] ] += thisBinValue;
	}

	//Delete the old distribution
	string name = LinearisedDistribution->GetName();
	string title = LinearisedDistribution->GetTitle();
	delete LinearisedDistribution;

	//Make the new distribution
	TH1F * delinearisedHistogram = new TH1F( name.c_str(), title.c_str(), binNumber - 2, distributionIndices->GetBinLowEdgesForRoot(0) );

	//Copy the data into the new distribution
	for ( unsigned int binIndex = 0; binIndex < binNumber; binIndex++ )
	{
		delinearisedHistogram->SetBinContent( binIndex, delinearisedDistribution[binIndex] );
	}

	//Return the new distribution
	return delinearisedHistogram;
}

vector< double > XvsYNormalisedFolding::DelineariseErrors( vector< double > LinearisedErrors )
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
vector< double > XvsYNormalisedFolding::CorrectedErrors()
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

//Return the names of the variables involved
vector< string > XvsYNormalisedFolding::VariableNames()
{
	vector< string > result( 1, xName );
	result.push_back( yName );
	return result;
}
