/**
  @interface IPlotMaker

  Gives all of the individual folding/unfolding plot classes a common interface, so that they call all be used by MonteCarloSummaryPlotMaker

  @author Benjamin M Wynne
  @date 08-03-2011
  */

#ifndef I_PLOTMAKER_H
#define I_PLOTMAKER_H

#include "InputNtuple.h"
#include "Distribution.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>

using namespace std;

class IPlotMaker
{
	public:
		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput ) = 0;
		virtual void StoreMiss( InputNtuple * TruthInput ) = 0;
		virtual void StoreFake( InputNtuple * ReconstructedInput ) = 0;
		virtual void StoreData( InputNtuple * DataInput ) = 0;

		//Do the unfolding
		virtual void Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing = false ) = 0;

		//Do a closure test
	        virtual bool ClosureTest( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing = false ) = 0;

		//Make a cross-check with MC
		virtual int MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing = false ) = 0;

		//Return a distribution for use in the cross-checks
		virtual Distribution * MonteCarloTruthForCrossCheck() = 0;

		//Return some plots
		virtual TH1F * CorrectedHistogram() = 0;
		virtual TH1F * UncorrectedHistogram() = 0;
		virtual TH1F * MCTruthHistogram() = 0;
		virtual TH2F * SmearingMatrix() = 0;

		//Copy the object
		virtual IPlotMaker * Clone( string NewPriorName ) = 0;

		//General info
		virtual string Description( bool WithSpaces ) = 0;
		virtual string PriorName() = 0;

		//Error info for corrected distribution
		virtual vector<double> CorrectedErrors() = 0;
};

#endif
