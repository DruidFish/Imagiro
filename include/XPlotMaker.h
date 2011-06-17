/**
  @class XPlotMaker

  Unfolds a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef X_PLOTMAKER_H
#define X_PLOTMAKER_H

#include "IUnfolder.h"
#include "IterativeUnfolding.h"
#include "IIndexCalculator.h"
#include <string>

using namespace std;

class XPlotMaker : public IUnfolder
{
	public:
		XPlotMaker();
		XPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum, double ScaleFactor = 1.0, bool Normalise = false );
		~XPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput );
		virtual void StoreMiss( IFileInput * TruthInput );
		virtual void StoreFake( IFileInput * ReconstructedInput );
		virtual void StoreData( IFileInput * DataInput );

		//Do the unfolding
		virtual void Unfold( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool SkipUnfolding = false, unsigned int ErrorMode = 0, bool WithSmoothing = false );

		//Do a closure test
		virtual bool ClosureTest( unsigned int MostIterations, double ChiSquaredThreshold, double KolmogorovThreshold, bool WithSmoothing = false );

		//Make a cross-check with MC
		virtual unsigned int MonteCarloCrossCheck( Distribution * ReferenceDistribution, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing = false );

		//Return a distribution for use in the cross-checks
		virtual Distribution * MonteCarloTruthForCrossCheck();

		//Return some plots
		virtual TH1F * CorrectedHistogram();
		virtual TH1F * UncorrectedHistogram();
		virtual TH1F * MCTruthHistogram();
		virtual TH2F * SmearingMatrix();

		//Copy the object
		virtual IUnfolder * Clone( string NewPriorName );

		//General info
		virtual string Description( bool WithSpaces );
		virtual string PriorName();

		//Error info for corrected distribution
		virtual vector< double > CorrectedErrors();
		virtual vector< double > DAgostiniErrors();
		virtual TH2F * DAgostiniCovariance();

		//Return the names of the variables involved
		virtual vector<string> VariableNames();

	private:
		//To be used with Clone
		XPlotMaker( string XVariableName, string PriorName, IIndexCalculator * DistributionIndices,
				unsigned int OriginalID, double ScaleFactor = 1.0, bool Normalise = false );

		unsigned int thisPlotID;
		IterativeUnfolding * XUnfolder;
		IIndexCalculator * distributionIndices;
		string xName, priorName;
		bool finalised, normalise;
		double scaleFactor;
		vector< double > correctedDataErrors, dagostiniErrors;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution;
		TH2F *smearingMatrix, *covarianceMatrix;
};

#endif
