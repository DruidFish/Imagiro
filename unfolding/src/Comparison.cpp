#include <sstream>
#include "Comparison.h"
#include "TH1F.h"
#include "TFile.h"
#include <iostream>

Comparison::Comparison() : internalID( 0 ), uniqueID( 0 ), name( "NoName" )
{
}

Comparison::Comparison( string Name, int UniqueID ) : internalID( 0 ), uniqueID( uniqueID ), name( Name )
{
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
	TH1F * firstPlot = FirstInput->MakeRootHistogram( firstPlotName, firstPlotName, false, Normalised );

	//Make a plot of the second distribution
	string secondPlotName = "SecondPlot" + internalName.str();
	TH1F * secondPlot = SecondInput->MakeRootHistogram( secondPlotName, secondPlotName, false, Normalised );

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
