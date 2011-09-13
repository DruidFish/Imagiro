/**
  @class MonteCarloSummaryPlotMaker

  Creates instances of a given "template" plot for each Monte Carlo sample, and then combines the results
  Systematic errors are given by the range of different resutls from using each Monte Carlo sample
  Cross-checks between samples are used to calculate convergence criteria for iterative unfolding

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef MONTE_CARLO_SUMMARY_PLOTMAKER_H
#define MONTE_CARLO_SUMMARY_PLOTMAKER_H

#include "MonteCarloInformation.h"
#include "IFileInput.h"
#include "IPlotMaker.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include <vector>
#include <string>

using namespace std;

class MonteCarloSummaryPlotMaker
{
	public:
		MonteCarloSummaryPlotMaker();
		MonteCarloSummaryPlotMaker( IPlotMaker * TemplatePlotMaker, MonteCarloInformation * PlotInformation, bool CombineMCMode = true );
		~MonteCarloSummaryPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput );
		void StoreMiss( IFileInput * TruthInput );
		void StoreFake( IFileInput * ReconstructedInput );
		void StoreData( IFileInput * DataInput );

		//Some plot formatting
		void SetYRange( double Minimum, double Maximum );
		void SetAxisLabels( string XAxis, string YAxis );
		void UseLogScale();

		//Do the calculation
		void Process( int ErrorMode = 0, bool WithSmoothing = false );

		//Save result to given file
		void SaveResult( TFile * OutputFile );

		//Return the names of the variables
		vector< string > VariableNames();

	private:
		int correctionType;
		vector< Distribution* > allTruthDistributions;
		vector< TH1F* > allTruthPlots;
		TH2F *smearingMatrix, *covarianceMatrix;
		MonteCarloInformation * mcInfo;
		TCanvas * plotCanvas;
		bool finalised, combineMode, manualRange, manualLabels, logScale;
		vector< IPlotMaker* > allPlots;
		vector< TH1F* > truthHistograms, reconstructedHistograms;
		TH1F *uncorrectedData, *correctedData, *statisticalErrors, *systematicErrors;
		double yRangeMinimum, yRangeMaximum;
		string dataDescription, xAxisLabel, yAxisLabel, plotDescription;
		vector< string > variableNames;
};

#endif
