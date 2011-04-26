/**
  @interface IFolder

  Gives all of the individual folding plot classes a common interface, so that they can all be used by MonteCarloSummaryPlotMaker

  @author Benjamin M Wynne
  @date 11-04-2011
  */

#ifndef I_FOLDER_H
#define I_FOLDER_H

#include "IPlotMaker.h"
#include "IFileInput.h"
#include "Distribution.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>

using namespace std;

class IFolder : public IPlotMaker
{
	public:
		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput ) = 0;
		virtual void StoreMiss( IFileInput * TruthInput ) = 0;
		virtual void StoreFake( IFileInput * ReconstructedInput ) = 0;
		virtual void StoreData( IFileInput * DataInput ) = 0;

		//Do the folding
		virtual void Fold() = 0;

		//Do a closure test
	        virtual bool ClosureTest() = 0;

		//Return some plots
		virtual TH1F * CorrectedHistogram() = 0;
		virtual TH1F * UncorrectedHistogram() = 0;
		virtual TH1F * MCTruthHistogram() = 0;
		virtual TH2F * SmearingMatrix() = 0;

		//Copy the object
		virtual IFolder * Clone( string NewPriorName ) = 0;

		//General info
		virtual string Description( bool WithSpaces ) = 0;
		virtual string PriorName() = 0;

		//Error info for corrected distribution
		virtual vector<double> CorrectedErrors() = 0;
};

#endif
