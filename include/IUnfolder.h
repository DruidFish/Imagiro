/**
  @interface IUnfolder

  Gives all of the individual unfolding plot classes a common interface, so that they can all be used by MonteCarloSummaryPlotMaker

  @author Benjamin M Wynne
  @date 08-03-2011
  */

#ifndef I_UNFOLDER_H
#define I_UNFOLDER_H

#include "IPlotMaker.h"
#include "IFileInput.h"
#include "Distribution.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>

using namespace std;

class IUnfolder : public IPlotMaker
{
	public:
		//Destructor
		virtual ~IUnfolder()
		{
		}

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput ) = 0;
		virtual void StoreMiss( IFileInput * TruthInput ) = 0;
		virtual void StoreFake( IFileInput * ReconstructedInput ) = 0;
		virtual void StoreData( IFileInput * DataInput ) = 0;

		//Do the unfolding
		virtual void Unfold( int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool SkipUnfolding = false, int ErrorMode = 0, bool WithSmoothing = false ) = 0;

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
		virtual IUnfolder * Clone( string NewPriorName ) = 0;

		//General info
		virtual string Description( bool WithSpaces ) = 0;
		virtual string PriorName() = 0;

		//Error info for corrected distribution
		virtual vector<double> CorrectedErrors() = 0;
		virtual vector<double> DAgostiniErrors() = 0;
		virtual TH2F * DAgostiniCovariance() = 0;
};

#endif
