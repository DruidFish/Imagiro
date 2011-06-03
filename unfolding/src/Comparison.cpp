/**
  @class Comparison

  An extremely simple class to return the Chi2 and Kolmogorov-Smirnoff comparison values betweeen two histograms
  The histogram naming scheme prevents the tedious complaints from Root about making multiple objects with the same name

  @author Benjamin M Wynne bwynne@cern.ch
  @date 08-03-2011
 */

#include <sstream>
#include "Comparison.h"
#include "TH1F.h"
#include "TFile.h"
#include <iostream>
#include <cmath>

const bool CLOSURE_FILE_OUTPUT = false;

Comparison::Comparison()
{
}

Comparison::Comparison( string Name, int UniqueID )
{
	internalID = 0;
	uniqueID = UniqueID;
	name = Name;
}

Comparison::~Comparison()
{
}

void Comparison::CompareDistributions( Distribution * FirstInput, Distribution * SecondInput, double & ChiSquared, double & Kolmogorov, bool Normalised, bool IsClosureTest )
{
	//Set up the name for both plots in this comparison
	stringstream internalName;
	internalName << name << "." << uniqueID << "." << internalID;

	TFile * debugFile;
	if ( IsClosureTest && CLOSURE_FILE_OUTPUT )
	{
		string fileName = internalName.str() + ".ClosureTest.root";
		debugFile = new TFile( fileName.c_str(), "RECREATE" );
	}

	//Make a plot of the first distribution
	string firstPlotName = "FirstPlot" + internalName.str();
	TH1F * firstPlot = FirstInput->MakeRootHistogram( firstPlotName, firstPlotName, Normalised );

	//Make a plot of the second distribution
	string secondPlotName = "SecondPlot" + internalName.str();
	TH1F * secondPlot = SecondInput->MakeRootHistogram( secondPlotName, secondPlotName, Normalised );

	if ( IsClosureTest && CLOSURE_FILE_OUTPUT )
	{
		firstPlot->Write();
		secondPlot->Write();
	}

	//Do the comparison
	ChiSquared = firstPlot->Chi2Test( secondPlot, "UUCHI2" );
	Kolmogorov = firstPlot->KolmogorovTest( secondPlot, "" );

	//For closure tests, give some more detailed info
	if ( IsClosureTest )
	{
		int maxDeviationIndex = 0;
		double sumErrors = 0;
		double maxDeviation = 0;
		double maxError = 0;
		for ( int binIndex = 1; binIndex <= firstPlot->GetNbinsX(); binIndex++ )
		{
			double error = secondPlot->GetBinContent( binIndex ) / firstPlot->GetBinContent( binIndex );
			if ( isnan( error ) )
			{
				error = 1.0;
			}
			if ( isinf( error ) )
			{
				error = 0.0;
			}
			sumErrors += error;

			double deviation = fabs( error - 1.0 );
			if ( deviation > maxDeviation || binIndex == 1 )
			{
				maxDeviation = deviation;
				maxError = error;
				maxDeviationIndex = binIndex;
			}
		}
		cout << "Average ratio ( unfolded bin ) / ( reference bin ) = " << sumErrors / (double)firstPlot->GetNbinsX() << endl;

		if ( maxDeviation > 1E-10 )
		{
			cout << "Greatest discrepancy ( " << maxError << " ) is in bin " << maxDeviationIndex << " of " << firstPlot->GetNbinsX() << endl;
		}
	}

	//Clean up
	if ( IsClosureTest && CLOSURE_FILE_OUTPUT )
	{
		debugFile->Close();
		delete debugFile;
	}
	else
	{
		delete firstPlot;
		delete secondPlot;
	}
	internalID++;
}
