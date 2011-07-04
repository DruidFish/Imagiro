/**
  @class UEdeltaPhiPlotMaker

  Unfolds a 1D distribution which takes multiple values within a single event
  UE-specific hard-coded addition to calculate delta-phi values

  @author Benjamin M Wynne bwynne@cern.ch
  @date 04-07-2011
 */

#ifndef UE_DELTA_PHI_PLOTMAKER_H
#define UE_DELTA_PHI_PLOTMAKER_H

#include "IUnfolder.h"
#include "IterativeUnfolding.h"
#include "IIndexCalculator.h"
#include <vector>
#include <string>

using namespace std;

class UEdeltaPhiPlotMaker : public IUnfolder
{
	public:
		UEdeltaPhiPlotMaker();
		UEdeltaPhiPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum,
				vector< string > OtherVariableNames = vector< string >(), double ScaleFactor = 1.0, bool Normalise = false );
		UEdeltaPhiPlotMaker( string XVariableName, string PriorName, vector< double > BinLowEdges,
				vector< string > OtherVariableNames = vector< string >(), double ScaleFactor = 1.0, bool Normalise = false );
		~UEdeltaPhiPlotMaker();

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
		UEdeltaPhiPlotMaker( vector< string > OtherVariableNames, string PriorName, IIndexCalculator * DistributionIndices,
				unsigned int OriginalID, double ScaleFactor, bool Normalise );

		//Connect up the sub event values
		void MakePairs( vector< double > * TruthValues, vector< double > * RecoValues );

		//Vectors to use in the pairing
		vector< unsigned int > unpairedTruth, unpairedReco;
		vector< pair< unsigned int, unsigned int > > truthRecoPairs;

		//UE-specific
		double PlusMinusPi( double NotInRange );

		unsigned int thisPlotID;
		IterativeUnfolding * XUnfolder;
		IIndexCalculator * distributionIndices;
		string xName, priorName;
		vector< string > otherPairingNames;
		bool finalised, normalise;
		double scaleFactor;
		vector< double > correctedDataErrors, dagostiniErrors;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution;
		TH2F *smearingMatrix, *covarianceMatrix;
};

#endif
