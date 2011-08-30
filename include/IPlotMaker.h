/**
  @interface IPlotMaker

  Gives all of the individual plot classes a common interface, so that they can all be used by MonteCarloSummaryPlotMaker

  @author Benjamin M Wynne
  @date 08-03-2011
  */

#ifndef I_PLOTMAKER_H
#define I_PLOTMAKER_H

#include "IFileInput.h"
#include "Distribution.h"
#include "SmearingMatrix.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>

using namespace std;

class IPlotMaker
{
	public:
		//Destructor
		virtual ~IPlotMaker()
		{
		}

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput ) = 0;
		virtual void StoreMiss( IFileInput * TruthInput ) = 0;
		virtual void StoreFake( IFileInput * ReconstructedInput ) = 0;
		virtual void StoreData( IFileInput * DataInput ) = 0;

		//Do the unfolding
		virtual void Correct( unsigned int MostIterations, bool SkipUnfolding = false, unsigned int ErrorMode = 0, bool WithSmoothing = false ) = 0;

		//Do a closure test
	        virtual bool ClosureTest( unsigned int MostIterations, bool WithSmoothing = false ) = 0;

		//Make a cross-check with MC
		virtual unsigned int MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing = false ) = 0;

		//Return a distribution for use in the cross-checks
		virtual Distribution * PriorDistributionForCrossCheck() = 0;
		virtual SmearingMatrix * SmearingMatrixForCrossCheck() = 0;

		//Return some plots
		virtual TH1F * CorrectedHistogram() = 0;
		virtual TH1F * UncorrectedHistogram() = 0;
		virtual TH1F * MCTruthHistogram() = 0;
		virtual TH2F * SmearingHistogram() = 0;

		//Copy the object
		virtual IPlotMaker * Clone( string NewPriorName ) = 0;

		//General info
		virtual string Description( bool WithSpaces ) = 0;
		virtual string PriorName() = 0;

		//Error info for corrected distribution
		virtual vector<double> CorrectedErrors() = 0;
		virtual TH2F * DAgostiniCovariance() = 0;

		//Return the names of the variables involved
		virtual vector<string> VariableNames() = 0;

		//Return the type of correction the plot will perform
		virtual int CorrectionMode() = 0;
};

#endif
