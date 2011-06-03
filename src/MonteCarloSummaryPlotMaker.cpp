/**
  @class MonteCarloSummaryPlotMaker

  Creates instances of a given "template" plot for each Monte Carlo sample, and then combines the results
  Systematic errors are given by the range of different resutls from using each Monte Carlo sample
  Cross-checks between samples are used to calculate convergence criteria for iterative unfolding

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#include "MonteCarloSummaryPlotMaker.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

const int MC_CHECK_OFFSET = 1;

//Default constructor - useless
MonteCarloSummaryPlotMaker::MonteCarloSummaryPlotMaker()
{
}

//Constructor for unfolding plots
MonteCarloSummaryPlotMaker::MonteCarloSummaryPlotMaker( IUnfolder * TemplatePlotMaker, MonteCarloInformation * PlotInformation, bool CombineMCMode )
{
	isUnfolding = true;
	finalised = false;
	manualRange = false;
	manualLabels = false;
	logScale = false;
	combineMode = CombineMCMode;
	mcInfo = PlotInformation;
	dataDescription = "";

	//Make a separate plot for each MC source
	bool usedTheTemplate = false;
	for ( unsigned int mcIndex = 0; mcIndex < mcInfo->NumberOfSources(); mcIndex++ )
	{
		string mcDescription = mcInfo->Description( mcIndex );

		//Don't make a redundant copy of the template
		if ( TemplatePlotMaker->PriorName() == mcDescription )
		{
			unfoldingPlots.push_back( TemplatePlotMaker );
			allPlots.push_back( TemplatePlotMaker );
			usedTheTemplate = true;
		}
		else
		{
			IUnfolder * clonePlot = TemplatePlotMaker->Clone( mcDescription );
			unfoldingPlots.push_back( clonePlot );
			allPlots.push_back( clonePlot );
		}

		//Make plots for testing the unfolding with MC
		crossCheckPlots.push_back( TemplatePlotMaker->Clone( mcDescription ) );
	}
	if ( !usedTheTemplate )
	{
		delete TemplatePlotMaker;
	}
}

//Constructor for folding plots
MonteCarloSummaryPlotMaker::MonteCarloSummaryPlotMaker( IFolder * TemplatePlotMaker, MonteCarloInformation * PlotInformation, bool CombineMCMode )
{
	isUnfolding = false;
	finalised = false;
	manualRange = false;
	manualLabels = false;
	logScale = false;
	combineMode = CombineMCMode;
	mcInfo = PlotInformation;
	dataDescription = "";

	//Make a separate plot for each MC source
	bool usedTheTemplate = false;
	for ( unsigned int mcIndex = 0; mcIndex < mcInfo->NumberOfSources(); mcIndex++ )
	{
		string mcDescription = mcInfo->Description( mcIndex );

		//Don't make a redundant copy of the template
		if ( TemplatePlotMaker->PriorName() == mcDescription )
		{
			foldingPlots.push_back( TemplatePlotMaker );
			allPlots.push_back( TemplatePlotMaker );
			usedTheTemplate = true;
		}
		else
		{
			IFolder * clonePlot = TemplatePlotMaker->Clone( mcDescription );
			foldingPlots.push_back( clonePlot );
			allPlots.push_back( clonePlot );
		}
	}
	if ( !usedTheTemplate )
	{
		delete TemplatePlotMaker;
	}
}

//Destructor
MonteCarloSummaryPlotMaker::~MonteCarloSummaryPlotMaker()
{
	if ( finalised )
	{
		for ( unsigned int truthIndex = 0; truthIndex < allTruthPlots.size(); truthIndex++ )
		{
			delete allTruthPlots[truthIndex];
		}
		delete smearingMatrix;
	}
	else
	{
		for ( unsigned int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			if ( isUnfolding )
			{
				delete crossCheckPlots[ plotIndex ];
			}
			delete foldingPlots[ plotIndex ];
		}
	}
}

//Take input values from ntuples
//To reduce file access, the appropriate row must already be in memory, the method does not change row
void MonteCarloSummaryPlotMaker::StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add matched MC events to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
	else
	{
		int inputIndex = TruthInput->DescriptionIndex();

		//Add the event to the smearing matrix if it's the correct MC input, or if MC inputs are combined
		if ( combineMode )
		{
			for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
			{
				allPlots[ mcIndex ]->StoreMatch( TruthInput, ReconstructedInput );
			}
		}
		else
		{
			allPlots[ inputIndex ]->StoreMatch( TruthInput, ReconstructedInput );
		}

		//Also store the event for the cross-check unfolding
		if ( isUnfolding )
		{
			if ( combineMode )
			{
				for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
				{
					crossCheckPlots[ mcIndex ]->StoreMatch( TruthInput, ReconstructedInput );
				}
			}
			else
			{
				crossCheckPlots[ inputIndex ]->StoreMatch( TruthInput, ReconstructedInput );
			}

			inputIndex += MC_CHECK_OFFSET;
			inputIndex %= allPlots.size();
			crossCheckPlots[ inputIndex ]->StoreData( ReconstructedInput );
		}
	}
}
void MonteCarloSummaryPlotMaker::StoreMiss( IFileInput * TruthInput )
{
	if ( finalised )
	{
		cerr << "Trying to add missed MC event to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
	else
	{
		int inputIndex = TruthInput->DescriptionIndex();

		//Add the event to the smearing matrix if it's the correct MC input, or if MC inputs are combined
		if ( combineMode )
		{
			for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
			{
				allPlots[ mcIndex ]->StoreMiss( TruthInput );
			}
		}       
		else
		{
			allPlots[ inputIndex ]->StoreMiss( TruthInput );
		}

		//Also store the event for the cross-check unfolding
		if ( isUnfolding )
		{
			if ( combineMode )
			{
				for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
				{
					crossCheckPlots[ mcIndex ]->StoreMiss( TruthInput );
				}
			}
			else
			{
				crossCheckPlots[ inputIndex ]->StoreMiss( TruthInput );
			}
		}
	}
}
void MonteCarloSummaryPlotMaker::StoreFake( IFileInput * ReconstructedInput )
{
	if ( finalised )
	{
		cerr << "Trying to add fake MC event to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		int inputIndex = ReconstructedInput->DescriptionIndex();

		//Add the event to the smearing matrix if it's the correct MC input, or if MC inputs are combined
		if ( combineMode )
		{
			for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
			{
				allPlots[ mcIndex ]->StoreFake( ReconstructedInput );
			}
		}
		else
		{
			allPlots[ inputIndex ]->StoreFake( ReconstructedInput );
		}

		//Also store the event for the cross-check unfolding
		if ( isUnfolding )
		{
			if ( combineMode )
			{
				for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
				{
					crossCheckPlots[ mcIndex ]->StoreFake( ReconstructedInput );
				}
			}
			else
			{
				crossCheckPlots[ inputIndex ]->StoreFake( ReconstructedInput );
			}

			inputIndex += MC_CHECK_OFFSET;
			inputIndex %= allPlots.size();
			crossCheckPlots[ inputIndex ]->StoreData( ReconstructedInput );
		}
	}
}
void MonteCarloSummaryPlotMaker::StoreData( IFileInput * DataInput )
{
	if ( finalised )
	{
		cerr << "Trying to add data event to finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}       
	else
	{
		for ( unsigned int mcIndex = 0; mcIndex < allPlots.size(); mcIndex++ )
		{
			allPlots[mcIndex]->StoreData( DataInput );
		}

		//Save the description of the data to use in setting the plot title
		if ( dataDescription.size() == 0 )
		{
			dataDescription = *( DataInput->Description() );
		}
	}
}

//Some plot formatting
void MonteCarloSummaryPlotMaker::SetYRange( double Minimum, double Maximum )
{
	if ( finalised )
	{
		cerr << "Trying to change plot range of finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
	else
	{
		if ( Minimum < Maximum )
		{
			yRangeMinimum = Minimum;
			yRangeMaximum = Maximum;
		}
		else
		{
			yRangeMinimum = Maximum;
			yRangeMaximum = Minimum;
		}

		manualRange = true;
	}
}
void MonteCarloSummaryPlotMaker::SetAxisLabels( string XAxis, string YAxis )
{
	if ( finalised )
	{
		cerr << "Trying to change axis labels of finalised MonteCarloSummaryPlotMaker" << endl;
		exit(1);
	}
	else
	{
		xAxisLabel = XAxis;
		yAxisLabel = YAxis;
		manualLabels = true;
	}
}
void MonteCarloSummaryPlotMaker::UseLogScale()
{
	logScale = true;
}

//Do the calculation
void MonteCarloSummaryPlotMaker::Process( int ErrorMode, bool WithSmoothing )
{
	if ( finalised )
	{
		cerr << "MonteCarloSummaryPlotMaker is already finalised" << endl;
		exit(1);
	}       
	else
	{
		int mostIterations = 0;
		double chiSquaredThreshold = 0.0;
		double kolmogorovThreshold = 0.0;
		string plotDescription = allPlots[0]->Description( true );

		//Do the unfolding cross-check to find out good conditions for convergence
		if ( isUnfolding )
		{
			cout << endl << "--------------- Started unfolding " << plotDescription << " ----------------" << endl;

			for ( unsigned int mcIndex = 0; mcIndex < crossCheckPlots.size(); mcIndex++ )
			{
				//Get the distribution that should be produced by the unfolding from MC truth
				int nextIndex = ( mcIndex + MC_CHECK_OFFSET ) % crossCheckPlots.size();
				Distribution * referenceDistribution = crossCheckPlots[ nextIndex ]->MonteCarloTruthForCrossCheck();

				//Unfold MC reco
				cout << endl << "Cross check - MC " << nextIndex << " reco with MC " << mcIndex << " prior" << endl;
				double convergenceChi2, convergenceKolmogorov;
				mostIterations += crossCheckPlots[ mcIndex ]->MonteCarloCrossCheck( referenceDistribution, convergenceChi2, convergenceKolmogorov );
				chiSquaredThreshold += convergenceChi2;
				kolmogorovThreshold += convergenceKolmogorov;
			}

			//Find the average values for the convergence criteria
			mostIterations = ceil( (double)mostIterations / (double)crossCheckPlots.size() );
			chiSquaredThreshold /= (double)crossCheckPlots.size();
			kolmogorovThreshold /= (double)crossCheckPlots.size();

			//Check that the iteration process is useful
			if ( mostIterations < 2 )
			{
				//Must apply the unfolding at least once, regardless
				mostIterations = 1;

				//Warn about potential problems
				cout << "Unfolding will only be applied once - no iteration. This suggests there is a problem: perhaps the smearing matrix is underpopulated?" << endl;
			}
			cout << "Chosen convergence criteria: max " << mostIterations << " iterations; chi2 below " << chiSquaredThreshold << "; KS above " << kolmogorovThreshold << endl;

			//Tidy up
			for ( unsigned int mcIndex = 0; mcIndex < crossCheckPlots.size(); mcIndex++ )
			{
				delete crossCheckPlots[ mcIndex ];
			}
		}
		else
		{
			cout << endl << "--------------- Started folding " << plotDescription << " ----------------" << endl;
		}

		//Perform closure tests
		int numberRemoved = 0;
		unsigned int numberFailed = 0;
		vector< bool > usePrior( allPlots.size(), true );
		for ( unsigned int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			cout << endl << "Closure test for " << mcInfo->Description( plotIndex ) << endl;
			bool closureWorked;
			if ( isUnfolding )
			{
				closureWorked = unfoldingPlots[plotIndex]->ClosureTest( mostIterations, chiSquaredThreshold, kolmogorovThreshold, WithSmoothing );
			}
			else
			{
				closureWorked = foldingPlots[plotIndex]->ClosureTest();
			}

			if ( !closureWorked )
			{
				numberFailed++;

				//if ( isUnfolding )
				//{
				//	//Remove priors that did not pass the closure test
				//	cout << "Removing " << mcInfo->Description( plotIndex ) << " from available priors" << endl;
				//	usePrior[ plotIndex ] = false;
				//	numberRemoved++;
				//}
			}
		}

		//Quit if too many prior distributions fail
		if ( numberFailed == allPlots.size() && isUnfolding )
		{
			cerr << "All priors failed their closure tests: something is really wrong here. Suggest you choose better binning / provide more MC stats" << endl;
			//exit(1);
		}
		else if ( numberFailed > (double)allPlots.size() / 2.0 )
		{
			cerr << "The majority of priors failed their closure tests. Suggest you choose better binning / provide more MC stats" << endl;
		}

		//Make a canvas to display the plots
		string plotName = allPlots[0]->Description(false) + "CorrectedDistribution";
		string plotTitle = allPlots[0]->Description(true) + " Corrected Distribution";
		plotCanvas = new TCanvas( plotName.c_str(), plotTitle.c_str(), 0, 0, 800, 600 );
		plotCanvas->Range( 0, 0, 1, 1 );
		plotCanvas->SetFillColor( kWhite );

		//Make a pad on the canvas for the main plot
		TPad * mainPad = new TPad( "mainPad", "mainPad", 0.01, 0.33, 0.99, 0.99 );
		mainPad->Draw();
		mainPad->cd();
		mainPad->SetTopMargin( 0.1 );
		mainPad->SetBottomMargin( 0.01 );
		mainPad->SetRightMargin( 0.1 );
		mainPad->SetFillStyle( 0 );
		mainPad->SetFillColor( kWhite );
		if ( logScale )
		{
			mainPad->SetLogy();
		}

		//Unfold each plot and retrieve the information
		vector<double> combinedCorrectedData, minimumCorrectedData, maximumCorrectedData, combinedStatisticErrors;
		TH1F *combinedCorrectedHistogramWithSystematics, *combinedCorrectedHistogramWithStatistics;
		allTruthPlots = vector< TH1F* >( allPlots.size(), NULL );
		bool firstPlot = true;
		for ( unsigned int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			//Unfold
			if ( isUnfolding )
			{
				if ( usePrior[ plotIndex ] )
				{
					cout << endl << "Unfolding " << allPlots[ plotIndex ]->Description(true) << " with " << allPlots[ plotIndex ]->PriorName() << endl;
				}
				else
				{
					cout << endl << "Skipping unfolding " << allPlots[ plotIndex ]->Description(true) << " with " << allPlots[ plotIndex ]->PriorName() << endl;
				}

				//Still need to run this to get the truth output
				unfoldingPlots[plotIndex]->Unfold( mostIterations, chiSquaredThreshold, kolmogorovThreshold, !usePrior[ plotIndex ], ErrorMode, WithSmoothing );
			}
			else
			{
				foldingPlots[ plotIndex ]->Fold();
			}

			//Make a local copy of the truth plot
			string truthPlotName = "localCopy" + allPlots[ plotIndex ]->PriorName() + "Truth";
			TH1F * newTruthHistogram = ( TH1F* )allPlots[ plotIndex ]->MCTruthHistogram()->Clone( truthPlotName.c_str() );
			allTruthPlots[ plotIndex ] = newTruthHistogram;

			//Only use the unfolded output if the closure test was passed
			if ( usePrior[ plotIndex ] )
			{
				//Get the error vector
				vector< double > plotErrors;
				if ( ErrorMode > 0 && isUnfolding )
				{
					plotErrors = unfoldingPlots[ plotIndex ]->DAgostiniErrors();

					vector< double > simpleErrors = allPlots[ plotIndex ]->CorrectedErrors();

					for ( unsigned int binIndex = 0; binIndex < plotErrors.size(); binIndex++ )
					{
						cout << plotErrors[ binIndex ] << ", " << simpleErrors[ binIndex ] << endl;
					}
				}
				else
				{
					plotErrors = allPlots[ plotIndex ]->CorrectedErrors();
				}

				//Load the corrected data into the combined distribution
				TH1F * correctedHistogram = allPlots[ plotIndex ]->CorrectedHistogram();
				for ( int binIndex = 0; binIndex < correctedHistogram->GetNbinsX() + 2; binIndex++ )
				{
					double binContent = correctedHistogram->GetBinContent( binIndex );

					if ( firstPlot )
					{
						//Initialise the distribution
						combinedCorrectedData.push_back( binContent );
						minimumCorrectedData.push_back( binContent );
						maximumCorrectedData.push_back( binContent );
						combinedStatisticErrors.push_back( plotErrors[ binIndex ] );
					}
					else
					{
						//Populate the distribution
						combinedCorrectedData[ binIndex ] += binContent;
						combinedStatisticErrors[ binIndex ] += plotErrors[ binIndex ];
						if ( binContent < minimumCorrectedData[ binIndex ] )
						{
							minimumCorrectedData[ binIndex ] = binContent;
						}
						if ( binContent > maximumCorrectedData[ binIndex ] )
						{
							maximumCorrectedData[ binIndex ] = binContent;
						}
					}
				}

				//Some things only need to be done once
				if ( firstPlot )
				{
					//Copy the format of the data histograms
					combinedCorrectedHistogramWithSystematics = new TH1F( *correctedHistogram );
					combinedCorrectedHistogramWithStatistics = new TH1F( *correctedHistogram );

					//Copy a smearing matrix if it exists
					smearingMatrix = NULL;
					if ( allPlots[plotIndex]->SmearingMatrix() )
					{
						string smearingName = allPlots[plotIndex]->Description(false) + "SmearingMatrix";
						string smearingTitle = allPlots[plotIndex]->Description(true) + " Smearing Matrix";
						smearingMatrix = ( TH2F* )allPlots[plotIndex]->SmearingMatrix()->Clone( smearingName.c_str() );
						smearingMatrix->SetTitle( smearingTitle.c_str() );
					}

					//Copy a covariance matrix
					if ( ErrorMode > 1 && isUnfolding )
					{
						string covarianceName = allPlots[plotIndex]->Description(false) + "CovarianceMatrix";
						string covarianceTitle = allPlots[plotIndex]->Description(true) + " Covariance Matrix";
						covarianceMatrix = ( TH2F* )unfoldingPlots[plotIndex]->DAgostiniCovariance()->Clone( covarianceName.c_str() );
						covarianceMatrix->SetTitle( covarianceTitle.c_str() );
					}
					else
					{
						covarianceMatrix = 0;
					}

					firstPlot = false;
				}
			}

			//Free some memory
			delete allPlots[ plotIndex ];
		}

		//Make a TGraph for the asymmetric errors
		int graphSize = combinedCorrectedData.size();
		Double_t xValues[ graphSize ];
		Double_t yValues[ graphSize ];
		Double_t xError[ graphSize ];
		Double_t yStatError[ graphSize ];
		Double_t yBothErrorLow[ graphSize ];
		Double_t yBothErrorHigh[ graphSize ];
		for ( int binIndex = 1; binIndex < graphSize - 1; binIndex++ )
		{
			//The x-axis values can just come straight from the histogram class
			xValues[ binIndex ] = combinedCorrectedHistogramWithStatistics->GetBinCenter( binIndex );
			xError[ binIndex ] = xValues[ binIndex ] - combinedCorrectedHistogramWithStatistics->GetBinLowEdge( binIndex );

			//Get the y-values by taking the means of the corrected distributions
			yValues[ binIndex ] = combinedCorrectedData[ binIndex ] / (double)( allPlots.size() - numberRemoved );
			combinedCorrectedHistogramWithStatistics->SetBinContent( binIndex, yValues[ binIndex ] );

			//Get the mean statistical error
			double statistic = combinedStatisticErrors[ binIndex ] / (double)( allPlots.size() - numberRemoved );
			yStatError[ binIndex ] = statistic;

			//The systematic errors come from the range of the corrected distributions
			double systematicLow = yValues[ binIndex ] - minimumCorrectedData[ binIndex ];
			double systematicHigh = maximumCorrectedData[ binIndex ] - yValues[ binIndex ];

			//Combine the statistical and systematic errors in quadrature
			yBothErrorLow[ binIndex ] = sqrt( ( systematicLow * systematicLow ) + ( statistic * statistic ) );
			yBothErrorHigh[ binIndex ] = sqrt( ( systematicHigh * systematicHigh ) + ( statistic * statistic ) );
		}

		//Get rid of the overflow bins
		double lowEdge = combinedCorrectedHistogramWithSystematics->GetBinLowEdge( 1 );
		double highEdge = combinedCorrectedHistogramWithSystematics->GetBinLowEdge( graphSize - 1 );
		double gap = highEdge - lowEdge;
		xValues[ 0 ] = lowEdge - gap;
		xValues[ graphSize - 1 ] = highEdge + gap;

		//Draw the combined error graph - asymmetric errors
		TGraphAsymmErrors * graphWithSystematics = new TGraphAsymmErrors( graphSize, xValues, yValues, xError, xError, yBothErrorLow, yBothErrorHigh );
		graphWithSystematics->SetName( "combinedErrorGraph" );
		graphWithSystematics->SetTitle( plotTitle.c_str() );
		graphWithSystematics->SetFillColor(30);
		graphWithSystematics->SetMarkerSize(1);
		graphWithSystematics->SetMarkerStyle(7);
		if ( manualRange )
		{
			//Manually set the y axis range
			graphWithSystematics->GetYaxis()->SetRangeUser( yRangeMinimum, yRangeMaximum );
		}
		graphWithSystematics->GetXaxis()->SetRangeUser( lowEdge, highEdge );
		graphWithSystematics->GetYaxis()->SetTitleOffset(0.95);
		if ( manualLabels )
		{
			//Maunally label the y axis
			graphWithSystematics->GetYaxis()->SetTitle( yAxisLabel.c_str() );
		}
		graphWithSystematics->Draw( "A2" );

		//Draw the stat error graph - symmetric errors
		TGraphErrors * graphWithStatistics = new TGraphErrors( graphSize, xValues, yValues, xError, yStatError );
		graphWithStatistics->SetName( "statErrorGraph" );
		graphWithStatistics->SetTitle( plotTitle.c_str() );
		graphWithStatistics->SetMarkerSize(0);
		graphWithStatistics->SetMarkerStyle(8);
		graphWithStatistics->Draw( "pz" );

		//Make the legend
		TLegend * lineColourKey = new TLegend( 0.55, 0.65, 0.95, 0.95 );
		lineColourKey->SetTextSize(0.04);
		lineColourKey->SetFillColor(kWhite);
		lineColourKey->SetBorderSize(0);
		lineColourKey->AddEntry( graphWithSystematics, dataDescription.c_str(), "lpf" );

		//Draw the MC truth histograms
		for ( unsigned int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			TH1F * truthPlot = allTruthPlots[ plotIndex ];
			truthPlot->SetLineColor( mcInfo->LineColour(plotIndex) );
			truthPlot->SetMarkerColor( mcInfo->LineColour(plotIndex) );
			truthPlot->SetLineStyle( mcInfo->LineStyle(plotIndex) );
			truthPlot->SetLineWidth(2.5);
			truthPlot->SetTitle( plotTitle.c_str() );
			truthPlot->Draw( "SAME" );

			//Add to legend
			string mcName = mcInfo->Description(plotIndex);
			if ( usePrior[ plotIndex ] )
			{
				mcName += " (prior)";
			}
			else
			{
				mcName += " (ref)";
			}
			lineColourKey->AddEntry( truthPlot, mcName.c_str(), "l" );
		}
		lineColourKey->Draw();

		//Make the ratio plot
		plotCanvas->cd();
		TPad * ratioPad = new TPad( "ratioPad", "ratioPad", 0.01, 0.01, 0.99, 0.32 );
		ratioPad->Draw();
		ratioPad->cd();
		ratioPad->SetTopMargin( 0.01 );
		ratioPad->SetBottomMargin( 0.3 );
		ratioPad->SetRightMargin( 0.1 );
		ratioPad->SetFillStyle( 0 );
		ratioPad->SetFillColor( kWhite );

		//All that is needed for the data is the combined error bars
		Double_t justOnes[ graphSize ];
		Double_t justZeros[ graphSize ];
		for ( int binIndex = 0; binIndex < graphSize; binIndex++ )
		{
			justOnes[ binIndex ] = 1.0;
			justZeros[ binIndex ] = 0.0;
			yBothErrorLow[ binIndex ] /= yValues[ binIndex ];
			yBothErrorHigh[ binIndex ] /= yValues[ binIndex ];
		}
		TGraphAsymmErrors * justDataErrors = new TGraphAsymmErrors( graphSize, xValues, justOnes, xError, xError, yBothErrorLow, yBothErrorHigh );
		justDataErrors->SetName( "ratioGraph" );
		justDataErrors->SetFillColor( 30 );
		justDataErrors->GetXaxis()->SetRangeUser( lowEdge, highEdge );
		justDataErrors->GetYaxis()->SetRangeUser( 0.75, 1.25 );
		justDataErrors->SetTitle();
		justDataErrors->GetYaxis()->SetTitleOffset( 0.5 );
		justDataErrors->GetYaxis()->SetTitle( "MC/Data" );
		justDataErrors->GetXaxis()->SetLabelSize( 0.11 );
		justDataErrors->GetYaxis()->SetLabelSize( 0.08 );
		justDataErrors->GetYaxis()->SetNdivisions( 408 );
		justDataErrors->GetXaxis()->SetTitleSize( 0.1 );
		justDataErrors->GetYaxis()->SetTitleSize( 0.1 );
		if ( manualLabels )
		{
			//Manually label the x axis
			justDataErrors->GetXaxis()->SetTitle( xAxisLabel.c_str() );
		}
		justDataErrors->Draw( "A2" );

		//Draw the stat error graph - symmetric errors
		TGraphErrors * justDataStatistics = new TGraphErrors( graphSize, xValues, justOnes, xError, justZeros );
		justDataStatistics->SetTitle( plotTitle.c_str() );
		justDataStatistics->SetMarkerStyle(1);
		justDataStatistics->Draw( "pz" );

		//Divide each truth plot by the data values
		for ( unsigned int plotIndex = 0; plotIndex < allPlots.size(); plotIndex++ )
		{
			//Make a copy of the existing plot
			string truthRatioPlotName = mcInfo->Description( plotIndex ) + "Ratio";
			TH1F * truthRatioPlot = ( TH1F* )allTruthPlots[ plotIndex ]->Clone( truthRatioPlotName.c_str() );

			//Divide by the data
			truthRatioPlot->Divide( combinedCorrectedHistogramWithStatistics );

			//Format
			truthRatioPlot->SetLineColor( mcInfo->LineColour( plotIndex ) );
			truthRatioPlot->SetMarkerColor( mcInfo->LineColour( plotIndex ) );
			truthRatioPlot->SetLineStyle( mcInfo->LineStyle( plotIndex ) );
			truthRatioPlot->SetLineWidth( 2.5 );
			truthRatioPlot->SetTitle( plotTitle.c_str() );
			truthRatioPlot->Draw( "SAME" );

			allTruthPlots.push_back( truthRatioPlot );
		}

		//Status message
		cout << endl << "--------------- Finished unfolding " << plotDescription << " ---------------" << endl;

		//Mark as done
		finalised = true;
	}
}

void MonteCarloSummaryPlotMaker::SaveResult( TFile * OutputFile )
{
	//Save the output
	OutputFile->cd();
	plotCanvas->Write();

	if ( smearingMatrix )
	{
		smearingMatrix->Write();
	}

	if ( covarianceMatrix )
	{
		covarianceMatrix->Write();
	}

	OutputFile->Flush();
}
