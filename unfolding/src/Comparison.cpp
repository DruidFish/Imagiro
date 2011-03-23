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
	if (IsClosureTest)
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

	if (IsClosureTest)
	{
		firstPlot->Write();
		secondPlot->Write();
	}

	//Do the comparison
	ChiSquared = firstPlot->Chi2Test( secondPlot, "UUCHI2" );
	Kolmogorov = firstPlot->KolmogorovTest( secondPlot, "" );

	//Clean up
	if (IsClosureTest)
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
