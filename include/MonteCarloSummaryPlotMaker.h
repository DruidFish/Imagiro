#ifndef MONTE_CARLO_SUMMARY_PLOTMAKER_H
#define MONTE_CARLO_SUMMARY_PLOTMAKER_H

#include "MonteCarloInformation.h"
#include "IPlotMaker.h"
#include "TCanvas.h"
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
		void Unfold( int MostIterations = 10, double ChiSquaredThreshold = 1.0, double KolmogorovThreshold = 1.0, bool WithSmoothing = false );

		//Return result
		TCanvas * ResultPlot();
		TH2F * SmearingMatrix();

	private:
		vector< TH1F* > allTruthPlots;
		TH2F * smearingMatrix;
		MonteCarloInformation * mcInfo;
		TCanvas * plotCanvas;
		bool finalised, combineMode;
		vector< IPlotMaker* > allPlots;
		double yRangeMinimum, yRangeMaximum;
};

#endif
