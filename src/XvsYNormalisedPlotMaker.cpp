/**
  @class XvsYNormalisedPlotMaker

  Unfolds a 2D distribution, and divides it by the unfolded 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#include "XvsYNormalisedPlotMaker.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

//Default constructor - useless
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker()
{
}

//Constructor with the names to use for the variables
XvsYNormalisedPlotMaker::XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
		int XBinNumber, double XMinimum, double XMaximum,
		int YBinNumber, double YMinimum, double YMaximum,
		double ScaleFactor ) : xName( XVariableName ), yName( YVariableName ), priorName( PriorName ), finalised( false ), scaleFactor(ScaleFactor)
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
	XUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + priorName, thisPlotID );

	//Store the y range
	minima.push_back( YMinimum );
	maxima.push_back( YMaximum );
	binNumbers.push_back( YBinNumber );

	//Make the x vs y unfolder
	XvsYUnfolder = new IterativeUnfolding( binNumbers, minima, maxima, xName + "vs" + yName + priorName, thisPlotID );

	//Set up the indices for the distributions
	DistributionIndices = new DataIndices( binNumbers, minima, maxima );

	//Set up the cross-check for data loss in delinearisation
	stringstream idString;
	idString << thisPlotID;
	string xvsyTruthName = xName + yName + priorName + "TruthCheck" + idString.str();
	string xTruthName = xName + priorName + "TruthCheck" + idString.str();
	xvsyTruthCheck = new TH1F( xvsyTruthName.c_str(), xvsyTruthName.c_str(), XBinNumber, XMinimum, XMaximum );
	xTruthCheck = new TH1F( xTruthName.c_str(), xTruthName.c_str(), XBinNumber, XMinimum, XMaximum );

	//Make a summary for the y data values
	yValueSummary = new StatisticsSummary();
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
IUnfolder * XvsYNormalisedPlotMaker::Clone( string NewPriorName )
{
	return new XvsYNormalisedPlotMaker( xName, yName, NewPriorName,
			DistributionIndices->GetBinNumber(0) - 2, DistributionIndices->GetMinima()[0], DistributionIndices->GetMaxima()[0],
			DistributionIndices->GetBinNumber(1) - 2, DistributionIndices->GetMinima()[1], DistributionIndices->GetMaxima()[1],
			scaleFactor );
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
void XvsYNormalisedPlotMaker::StoreMiss( IFileInput * TruthInput )
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
void XvsYNormalisedPlotMaker::StoreFake( IFileInput * ReconstructedInput )
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
void XvsYNormalisedPlotMaker::StoreData( IFileInput * DataInput )
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
		yValueSummary->StoreEvent( yDataValue, dataWeight );

		//Store values for performing the delinearisation
		DistributionIndices->StoreDataValue( dataValues, dataWeight );
	} 
}

//Do the unfolding
void XvsYNormalisedPlotMaker::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool SkipUnfolding, int ErrorMode, bool WithSmoothing )
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
			XUnfolder->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, ErrorMode, WithSmoothing );
			XvsYUnfolder->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, ErrorMode, WithSmoothing );
		}

		//Make some plot titles
		stringstream uniqueIDString;
		uniqueIDString << thisPlotID;
		string XFullName = xName + priorName + uniqueIDString.str();
		string XFullTitle = xName + " using " + priorName;
		string XvsYName = xName + "vs" + yName + priorName + uniqueIDString.str();
		string XvsYTitle = xName + " vs " + yName + " using " + priorName;

		//Retrieve the results
		TH1F * XCorrected = XUnfolder->GetUnfoldedHistogram( XFullName + "Corrected", XFullTitle + " Corrected Distribution" );
		TH1F * XvsYCorrected = XvsYUnfolder->GetUnfoldedHistogram( XvsYName + "Corrected", XvsYTitle + " Corrected Distribution" );

		//Retrieve some other bits for debug
		TH1F * XUncorrected = XUnfolder->GetUncorrectedDataHistogram( XFullName + "Uncorrected", XFullTitle + " Uncorrected Distribution" );
		TH1F * XvsYUncorrected = XvsYUnfolder->GetUncorrectedDataHistogram( XvsYName + "Uncorrected", XvsYTitle + " Uncorrected Distribution" );

		TH1F * XTruth = XUnfolder->GetTruthHistogram( XFullName + "Truth", XFullTitle + " Truth Distribution" );
		TH1F * XvsYTruth = XvsYUnfolder->GetTruthHistogram( XvsYName + "Truth", XvsYTitle + " Truth Distribution" );

		//TH2F * XSmearing = XUnfolder->GetSmearingMatrix( XFullName + "Smearing", XFullTitle + " Smearing Matrix" );
		//TH2F * XvsYSmearing = XvsYUnfolder->GetSmearingMatrix( XvsYName + "Smearing", XvsYTitle + " Smearing Matrix" );
		if ( ErrorMode > 1 )
		{
			covarianceMatrix = XvsYUnfolder->DAgostiniCovariance( XvsYName + "Covariance", XvsYTitle + " Covariance Matrix" );
		}

		//De-linearise the x vs y distributions
		TH1F * DelinearisedXvsYCorrected = Delinearise(XvsYCorrected);
		TH1F * DelinearisedXvsYUncorrected = Delinearise(XvsYUncorrected);
		TH1F * DelinearisedXvsYTruth = Delinearise(XvsYTruth);

		//Get the error vectors
		vector<double> XErrors = XUnfolder->SumOfDataWeightSquares();
		vector<double> XvsYErrors = XvsYUnfolder->SumOfDataWeightSquares();
		vector<double> XVariance = XUnfolder->DAgostiniVariance();
		vector<double> XvsYVariance = XvsYUnfolder->DAgostiniVariance();

		//Delinearise the x vs y errors
		vector<double> delinearisedXvsYErrors = DelineariseErrors( XvsYErrors );

		//Combine errors
		for ( unsigned int binIndex = 0; binIndex < XErrors.size(); binIndex++ )
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

		//Same for the D'Agostini variance
		if ( ErrorMode > 0 )
		{
			vector<double> delinearisedXvsYVariance = DelineariseErrors( XvsYVariance );
			for ( unsigned int binIndex = 0; binIndex < XVariance.size(); binIndex++ )
			{
				//Add the errors from the x vs y distribution and the divisor
				//The formula is stright from ROOT::TH1::Divide, so I hope it's right
				double XBinValue = XCorrected->GetBinContent(binIndex);
				double XvsYBinValue = DelinearisedXvsYCorrected->GetBinContent(binIndex);
				double componentOne = XVariance[binIndex] * XvsYBinValue * XvsYBinValue;
				double componentTwo = delinearisedXvsYVariance[binIndex] * XBinValue * XBinValue;
				double componentThree = XBinValue * XBinValue;
				double combinedError = ( sqrt( componentOne + componentTwo ) / componentThree );
				combinedError *= scaleFactor;
				dagostiniErrors.push_back( combinedError );
			}
		}

		//Normalise the x vs y distributions and scale appropriately
		DelinearisedXvsYCorrected->Divide( DelinearisedXvsYCorrected, XCorrected, scaleFactor, 1.0 );
		DelinearisedXvsYUncorrected->Divide( DelinearisedXvsYUncorrected, XUncorrected, scaleFactor, 1.0 );
		DelinearisedXvsYTruth->Divide( DelinearisedXvsYTruth, XTruth, scaleFactor, 1.0 );

		//Check for data loss in the delinearisation
		xvsyTruthCheck->Divide( xvsyTruthCheck, xTruthCheck, scaleFactor, 1.0 );
		bool dataLost = false;
		double averagePercentError = 0.0;
		for ( int binIndex = 0; binIndex < xTruthCheck->GetNbinsX() + 2; binIndex++ )
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
				if ( ErrorMode > 0 )
				{
					dagostiniErrors[ binIndex ] *= ( 1.0 + errorFraction );
				}

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
		averagePercentError /= (double)( xTruthCheck->GetNbinsX() + 2 );
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
		delete XCorrected;
		delete XUncorrected;
		delete XTruth;
		delete yValueSummary;

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
		smearingMatrix = new TH2F( XvsYName.c_str(), XvsYTitle.c_str(), 1, 0.0, 1.0, 1, 0.0, 1.0 );//XvsYSmearing;
		string smearingXTitle = xName + " vs " + yName + " Truth";
		string smearingYTitle = xName + " vs " + yName + " Reconstructed";
		smearingMatrix->SetXTitle( smearingXTitle.c_str() );
		smearingMatrix->SetYTitle( smearingYTitle.c_str() );

		if ( ErrorMode > 1 )
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
bool XvsYNormalisedPlotMaker::ClosureTest( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	bool result = XvsYUnfolder->ClosureTest( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
	result &= XUnfolder->ClosureTest( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
	return result;
}

//Make a cross-check with MC
int XvsYNormalisedPlotMaker::MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing )
{
	return XvsYUnfolder->MonteCarloCrossCheck( ReferenceDistribution, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );
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

//Return a distribution for use in the cross-checks
Distribution * XvsYNormalisedPlotMaker::MonteCarloTruthForCrossCheck()
{
	return XvsYUnfolder->GetTruthDistribution();
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
vector<double> XvsYNormalisedPlotMaker::DAgostiniErrors()
{
	if ( finalised )
	{
		return dagostiniErrors;
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
