#include "MonteCarloSummaryPlotMaker.h"
#include "TLegend.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

//Default constructor - useless
MonteCarloSummaryPlotMaker::MonteCarloSummaryPlotMaker()
{
}

//Constructor with the names to use for the variables
MonteCarloSummaryPlotMaker::MonteCarloSummaryPlotMaker( IPlotMaker * TemplatePlotMaker, MonteCarloInformation * PlotInformation,
		double YMinimum, double YMaximum, bool CombineMCMode ) : finalised( false ), combineMode(CombineMCMode), mcInfo(PlotInformation), yRangeMinimum( YMinimum ), yRangeMaximum( YMaximum )
{
	//Make a separate plot for each MC source
	for ( int mcIndex = 0; mcIndex < mcInfo->NumberOfSources(); mcIndex++ )
	{
		string mcDescription = mcInfo->Description( mcIndex );
		if ( TemplatePlotMaker->PriorName() == mcDescription )
		{
			allPlots.push_back( TemplatePlotMaker );
		}
		else
		{
			IPlotMaker * thisMCPlot = TemplatePlotMaker->Clone( mcDescription );
			allPlots.push_back( thisMCPlot );
		}
	}
}

//Destructor
MonteCarloSummaryPlotMaker::~MonteCarloSummaryPlotMaker()
{
	if ( finalised )
	{
		allTruthPlots.clear();
		delete smearingMatrix;
	}
	else
	{
		allPlots.clear();
	}
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void MonteCarloSummaryPlotMaker::StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
	else
	{
		for ( int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
		{
			if ( combineMode || mcInfo->Description( mcIndex ) == *( TruthInput->Description() ) )
			{
				allPlots[mcIndex]->StoreMatch( TruthInput, ReconstructedInput );
			}
		}
	}
}
void MonteCarloSummaryPlotMaker::StoreMiss( InputNtuple * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		for ( int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
		{
			if ( combineMode || mcInfo->Description( mcIndex ) == *( TruthInput->Description() ) )
			{
				allPlots[mcIndex]->StoreMiss( TruthInput );
			}
		}
	}
}
void MonteCarloSummaryPlotMaker::StoreFake( InputNtuple * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		for ( int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
		{
			if ( combineMode || mcInfo->Description( mcIndex ) == *( ReconstructedInput->Description() ) )
			{
				allPlots[mcIndex]->StoreFake( ReconstructedInput );
			}
		}
	}
}
void MonteCarloSummaryPlotMaker::StoreData( InputNtuple * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		for ( int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
		{
			allPlots[mcIndex]->StoreData( DataInput );
		}
	}
}

//Do the unfolding
void MonteCarloSummaryPlotMaker::Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "MonteCarloSummaryPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		//Make a canvas to display the plots
		string plotName = allPlots[0]->Description(false) + "CorrectedDistribution";
		string plotTitle = allPlots[0]->Description(true) + " Corrected Distribution";
		plotCanvas = new TCanvas( plotName.c_str(), plotTitle.c_str() );
		plotCanvas->SetFillColor(kWhite);

		//Unfold each plot and retrieve the information
		vector<double> combinedCorrectedData, minimumCorrectedData, maximumCorrectedData, combinedStatisticErrors;
		TH1F *combinedCorrectedHistogramWithSystematics, *combinedCorrectedHistogramWithStatistics;
		for ( int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			//Unfold
			cout << endl << "Unfolding " << allPlots[plotIndex]->Description(true) << " with " << allPlots[plotIndex]->PriorName() << endl;
			allPlots[plotIndex]->Unfold( MostIterations, ChiSquaredThreshold, KolmogorovThreshold, WithSmoothing );

			//Get the error vector
			vector<double> plotErrors = allPlots[ plotIndex ]->CorrectedErrors();

			//Make a local copy of the truth plot
			string truthPlotName = "localCopy" + allPlots[plotIndex]->PriorName() + "Truth";
			TH1F * newTruthPlot = ( TH1F* )allPlots[ plotIndex ]->MCTruthDistribution()->Clone( truthPlotName.c_str() );
			allTruthPlots.push_back(newTruthPlot);

			//Load the corrected data into the combined distribution
			TH1F * correctedDistribution = allPlots[ plotIndex ]->CorrectedDistribution();
			for ( int binIndex = 0; binIndex < correctedDistribution->GetNbinsX() + 2; binIndex++ )
			{
				double binContent = correctedDistribution->GetBinContent( binIndex );

				if ( plotIndex == 0 )
				{
					//Initialise the distribution
					combinedCorrectedData.push_back( binContent );
					minimumCorrectedData.push_back( binContent );
					maximumCorrectedData.push_back( binContent );
					combinedStatisticErrors.push_back( plotErrors[binIndex] );
				}
				else
				{
					//Populate the distribution
					combinedCorrectedData[binIndex] += binContent;
					combinedStatisticErrors[binIndex] += plotErrors[binIndex];
					if ( binContent < minimumCorrectedData[binIndex] )
					{
						minimumCorrectedData[binIndex] = binContent;
					}
					if ( binContent > maximumCorrectedData[binIndex] )
					{
						maximumCorrectedData[binIndex] = binContent;
					}
				}
			}

			//Some things only need to do once
			if ( plotIndex == 0 )
			{
				//Copy the format of the data histograms
				combinedCorrectedHistogramWithSystematics = new TH1F( *correctedDistribution );
				combinedCorrectedHistogramWithStatistics = new TH1F( *correctedDistribution );

				//Copy a smearing matrix
				string smearingName = allPlots[plotIndex]->Description(false) + "SmearingMatrix";
				string smearingTitle = allPlots[plotIndex]->Description(true) + " Smearing Matrix";
				smearingMatrix = ( TH2F* )allPlots[plotIndex]->SmearingMatrix()->Clone( smearingName.c_str() );
				smearingMatrix->SetTitle( smearingTitle.c_str() );
				smearingMatrix->SetStats(false);
			}

			//Free some memory
			delete allPlots[ plotIndex ];
		}

		//Make histograms of the combined corrected distributions
		for ( int binIndex = 0; binIndex < combinedCorrectedData.size(); binIndex++ )
		{
			//Take the mean of the central values
			combinedCorrectedData[binIndex] /= allPlots.size();

			//The systematic error is just the range of the corrected data
			combinedCorrectedHistogramWithSystematics->SetBinContent( binIndex, combinedCorrectedData[binIndex] );
			double systematic = ( maximumCorrectedData[binIndex] - minimumCorrectedData[binIndex] ) / 2.0;
			double statistic = combinedStatisticErrors[binIndex] / (double)allPlots.size();
			double systematicAndStatisticError = sqrt( ( systematic * systematic ) + ( statistic * statistic ) );
			combinedCorrectedHistogramWithSystematics->SetBinError( binIndex, systematicAndStatisticError );

			//Get the statistical error from one of the plotmakers
			combinedCorrectedHistogramWithStatistics->SetBinContent( binIndex, combinedCorrectedData[binIndex] );
			combinedCorrectedHistogramWithStatistics->SetBinError( binIndex, statistic );
		}

		//Draw the data histogram with systematic errors
		combinedCorrectedHistogramWithSystematics->SetTitle( plotTitle.c_str() );
		combinedCorrectedHistogramWithSystematics->SetStats(false);
		combinedCorrectedHistogramWithSystematics->SetMarkerStyle(7);
		combinedCorrectedHistogramWithSystematics->SetMarkerSize(1);
		combinedCorrectedHistogramWithSystematics->SetFillColor(30);
		if ( yRangeMinimum != yRangeMaximum )
		{
			//Manually set the y axis range
			if ( yRangeMinimum < yRangeMaximum )
			{
				combinedCorrectedHistogramWithSystematics->GetYaxis()->SetRangeUser( yRangeMinimum, yRangeMaximum );
			}
			else
			{
				combinedCorrectedHistogramWithSystematics->GetYaxis()->SetRangeUser( yRangeMaximum, yRangeMinimum );
			}
		}
		combinedCorrectedHistogramWithSystematics->Draw( "E2" );

		//Draw the data histogram with statistical errors
		combinedCorrectedHistogramWithStatistics->SetTitle( plotTitle.c_str() );
		combinedCorrectedHistogramWithStatistics->SetStats(false);
		combinedCorrectedHistogramWithStatistics->SetMarkerStyle(8);
		combinedCorrectedHistogramWithStatistics->SetMarkerSize(0);
		combinedCorrectedHistogramWithStatistics->Draw( "SAME" );

		//Make the legend
		TLegend * lineColourKey = new TLegend( 0.25, 0.05, 0.50, 0.23 );
		lineColourKey->SetTextSize(0.04);
		lineColourKey->SetFillColor(kWhite);
		lineColourKey->SetBorderSize(0);
		lineColourKey->AddEntry( combinedCorrectedHistogramWithSystematics, "Data 2010", "lpf" );

		//Draw the MC truth histograms
		for ( int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			TH1F * truthPlot = allTruthPlots[ plotIndex ];
			truthPlot->SetLineColor( mcInfo->LineColour(plotIndex) );
			truthPlot->SetMarkerColor( mcInfo->LineColour(plotIndex) );
			truthPlot->SetLineStyle( mcInfo->LineStyle(plotIndex) );
			truthPlot->SetLineWidth(2.5);
			truthPlot->SetTitle( plotTitle.c_str() );
			truthPlot->SetStats(false);
			truthPlot->Draw( "SAME" );

			//Add to legend
			lineColourKey->AddEntry( truthPlot, mcInfo->Description(plotIndex).c_str(), "l" );
		}
		lineColourKey->Draw();

		//Mark as done
		finalised = true;
	}
}

//Return result
TCanvas * MonteCarloSummaryPlotMaker::ResultPlot()
{
	if ( finalised )
	{
		return plotCanvas;
	}
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
}
TH2F * MonteCarloSummaryPlotMaker::SmearingMatrix()
{
	if ( finalised )
	{
		return smearingMatrix;
	}       
	else
	{
		cerr << "Trying to retrieve smearing matrix from unfinalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
}
