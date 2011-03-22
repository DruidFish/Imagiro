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
#include "IPlotMaker.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <vector>
#include <string>

using namespace std;

class MonteCarloSummaryPlotMaker
{
	public:
		MonteCarloSummaryPlotMaker();
		MonteCarloSummaryPlotMaker( IPlotMaker * TemplatePlotMaker, MonteCarloInformation * PlotInformation, double YMinimum = 0.0, double YMaximum = 0.0, bool CombineMCMode = true );
		~MonteCarloSummaryPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		void StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput );
		void StoreMiss( InputNtuple * TruthInput );
		void StoreFake( InputNtuple * ReconstructedInput );
		void StoreData( InputNtuple * DataInput );

		//Do the unfolding
		void Unfold( bool WithSmoothing = false );

		//Return result
		TCanvas * ResultPlot();
		TH2F * SmearingMatrix();

		//Create an ATLAS style object
		TStyle * AtlasStyle( string Name );

	private:
		vector< Distribution* > allTruthDistributions;
		vector< TH1F* > allTruthPlots;
		TH2F * smearingMatrix;
		MonteCarloInformation * mcInfo;
		TCanvas * plotCanvas;
		bool finalised, combineMode;
		vector< IPlotMaker* > allPlots;
		vector< IPlotMaker* > crossCheckPlots;
		double yRangeMinimum, yRangeMaximum;
		string dataDescription;
};

#endif
