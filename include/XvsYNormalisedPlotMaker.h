/**
  @class XvsYNormalisedPlotMaker

  Unfolds a 2D distribution, and divides it by the normalised 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#ifndef X_VS_Y_NORMALISED_PLOTMAKER_H
#define X_VS_Y_NORMALISED_PLOTMAKER_H

#include "IPlotMaker.h"
#include "StatisticsSummary.h"
#include "IterativeUnfolding.h"
#include "DataIndices.h"
#include <string>

using namespace std;

class XvsYNormalisedPlotMaker : public IPlotMaker
{
	public:
		XvsYNormalisedPlotMaker();
		XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName, int XBinNumber, double XMinimum, double XMaximum,
				int YBinNumber, double YMinimum, double YMaximum, double ScaleFactor = 1.0 );
		~XvsYNormalisedPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput );
		virtual void StoreMiss( InputNtuple * TruthInput );
		virtual void StoreFake( InputNtuple * ReconstructedInput );
		virtual void StoreData( InputNtuple * DataInput );

		//Do the unfolding
		virtual void Unfold( int MostIterations = 10, double ChiSquaredThreshold = 1.0, double KolmogorovThreshold = 1.0, bool WithSmoothing = false );

		//Make a cross-check with MC
		virtual int MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing = false );

		//Return a distribution for use in the cross-checks
		virtual Distribution * MonteCarloTruthForCrossCheck();

		//Return some plots
		virtual TH1F * CorrectedHistogram();
		virtual TH1F * UncorrectedHistogram();
		virtual TH1F * MCTruthHistogram();
		virtual TH2F * SmearingMatrix();

		//Copy the object
		virtual IPlotMaker * Clone( string NewPriorName );

		//General info
		virtual string Description( bool WithSpaces );
		virtual string PriorName();

		//Error info for corrected distribution
		virtual vector<double> CorrectedErrors();

	private:
		//WARNING: this method deletes the argument object
		TH1F * Delinearise( TH1F * LinearisedDistribution );
		vector<double> DelineariseErrors( vector<double> InputSumWeightSquares );

		int thisPlotID;
		IterativeUnfolding *XvsYUnfolder, *XUnfolder;
		DataIndices *DistributionIndices;
		string xName, yName, priorName;
		bool finalised;
		double scaleFactor;
		vector<double> correctedDataErrors;
		StatisticsSummary * yValueSummary;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution, *xvsyTruthCheck, *xTruthCheck;
		TH2F *smearingMatrix;
};

#endif
