/**
  @class XFolding

  Applies the smearing matrix to a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef X_FOLDING_H
#define X_FOLDING_H

#include "IPlotMaker.h"
#include "Folding.h"
#include "Indices.h"
#include <string>

using namespace std;

class XFolding : public IPlotMaker
{
	public:
		XFolding();
		XFolding( string XVariableName, string PriorName, int XBinNumber, double XMinimum, double XMaximum, double ScaleFactor = 1.0 );
		~XFolding();

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
		int thisPlotID;
		Folding * XFolder;
		Indices * DistributionIndices;
		string xName, priorName;
		bool finalised;
		double scaleFactor;
		vector<double> correctedInputErrors;
		TH1F *foldedDistribution, *inputDistribution, *reconstructedDistribution;
		TH2F *smearingMatrix;
};

#endif
