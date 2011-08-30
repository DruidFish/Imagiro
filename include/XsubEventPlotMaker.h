/**
  @class XsubEventPlotMaker

  Unfolds a 1D distribution which takes multiple values within a single event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 03-07-2011
 */

#ifndef X_SUB_EVENT_PLOTMAKER_H
#define X_SUB_EVENT_PLOTMAKER_H

#include "IPlotMaker.h"
#include "ICorrection.h"
#include "IIndexCalculator.h"

using namespace std;

class XsubEventPlotMaker : public IPlotMaker
{
	public:
		XsubEventPlotMaker();
		XsubEventPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum,
				int CorrectionMode = 2, vector< string > OtherVariableNames = vector< string >(), double ScaleFactor = 1.0, bool Normalise = false );
		XsubEventPlotMaker( string XVariableName, string PriorName, vector< double > BinLowEdges,
				int CorrectionMode = 2, vector< string > OtherVariableNames = vector< string >(), double ScaleFactor = 1.0, bool Normalise = false );

		virtual ~XsubEventPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput );
		virtual void StoreMiss( IFileInput * TruthInput );
		virtual void StoreFake( IFileInput * ReconstructedInput );
		virtual void StoreData( IFileInput * DataInput );

		//Do the unfolding
		virtual void Correct( unsigned int MostIterations, bool SkipUnfolding = false, unsigned int ErrorMode = 0, bool WithSmoothing = false );

		//Do a closure test
		virtual bool ClosureTest( unsigned int MostIterations, bool WithSmoothing = false );

		//Make a cross-check with MC
		virtual unsigned int MonteCarloCrossCheck( Distribution * InputPriorDistribution, SmearingMatrix * InputSmearing, bool WithSmoothing = false );

		//Return a distribution for use in the cross-checks
		virtual Distribution * PriorDistributionForCrossCheck();
		virtual SmearingMatrix * SmearingMatrixForCrossCheck();

		//Return some plots
		virtual TH1F * CorrectedHistogram();
		virtual TH1F * UncorrectedHistogram();
		virtual TH1F * MCTruthHistogram();
		virtual TH2F * SmearingHistogram();

		//Copy the object
		virtual XsubEventPlotMaker * Clone( string NewPriorName );

		//General info
		virtual string Description( bool WithSpaces );
		virtual string PriorName();

		//Error info for corrected distribution
		virtual vector< double > CorrectedErrors();
		virtual TH2F * DAgostiniCovariance();

		//Return the names of the variables involved
		virtual vector< string > VariableNames();

		//Return the type of correction the plot will perform
		virtual int CorrectionMode();

	private:
		//To be used with Clone
		XsubEventPlotMaker( vector< string > OtherVariableNames, string PriorName, IIndexCalculator * DistributionIndices,
				unsigned int OriginalID, int CorrectionMode, double ScaleFactor, bool Normalise );

		//Connect up the sub event values
		void MakePairs( vector< double > * TruthValues, vector< double > * RecoValues );

		//Vectors to use in the pairing
		vector< unsigned int > unpairedTruth, unpairedReco;
		vector< pair< unsigned int, unsigned int > > truthRecoPairs;

		int correctionType;
		unsigned int thisPlotID;
		ICorrection * XUnfolder;
		IIndexCalculator * distributionIndices;
		string xName, priorName;
		vector< string > otherPairingNames;
		bool finalised, normalise;
		double scaleFactor;
		vector< double > correctedDataErrors;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution;
		TH2F *smearingMatrix, *covarianceMatrix;
};

#endif
