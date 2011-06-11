/**
  @interface IPlotMaker

  Gives all of the individual folding/unfolding plot classes a common interface, so that they can all be used by MonteCarloSummaryPlotMaker

  @author Benjamin M Wynne
  @date 08-03-2011
  */

#ifndef I_PLOTMAKER_H
#define I_PLOTMAKER_H

#include "IFileInput.h"
#include "Distribution.h"
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

		//Return some plots
		virtual TH1F * CorrectedHistogram() = 0;
		virtual TH1F * UncorrectedHistogram() = 0;
		virtual TH1F * MCTruthHistogram() = 0;
		virtual TH2F * SmearingMatrix() = 0;

		//General info
		virtual string Description( bool WithSpaces ) = 0;
		virtual string PriorName() = 0;

		//Error info for corrected distribution
		virtual vector<double> CorrectedErrors() = 0;

		//Return the names of the variables involved
		virtual vector<string> VariableNames() = 0;
};

#endif
