/**
  @class XPlotMaker

  Unfolds a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef X_PLOTMAKER_H
#define X_PLOTMAKER_H

#include "IPlotMaker.h"
#include "ICorrection.h"
#include "IIndexCalculator.h"
#include "TRandom3.h"
#include <string>

using namespace std;

class XPlotMaker : public IPlotMaker
{
	public:
		XPlotMaker();
		XPlotMaker( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum, int CorrectionMode = 2, double ScaleFactor = 1.0, bool Normalise = false );
		XPlotMaker( string XVariableName, string PriorName, vector< double > BinLowEdges, int CorrectionMode = 2, double ScaleFactor = 1.0, bool Normalise = false );

		virtual ~XPlotMaker();

		//Set up a systematic error study
		virtual void AddSystematic( vector< double > SystematicOffset, vector< double > SystematicWidth, unsigned int NumberOfPseudoExperiments );

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
		virtual TH1F * MCRecoHistogram();
		virtual TH2F * SmearingHistogram();
		virtual vector< TH1F* > SystematicHistograms();

		//Copy the object
		virtual XPlotMaker * Clone( string NewPriorName );

		//General info
		virtual string Description( bool WithSpaces );
		virtual string PriorName();

		//Error info for corrected distribution
		virtual vector< double > CorrectedErrors();
		virtual TH2F * DAgostiniCovariance();

		//Return the names of the variables involved
		virtual vector<string> VariableNames();

		//Return the type of correction the plot will perform
                virtual int CorrectionMode();

	private:
		//To be used with Clone
		XPlotMaker( string XVariableName, string PriorName, IIndexCalculator * DistributionIndices,
				unsigned int OriginalID, int CorrectionMode, double ScaleFactor, bool Normalise, vector<double> InputOffsets, vector<double> InputWidths );

		//Instantiate the corrector
		ICorrection * MakeCorrector( int CorrectionMode );

		int correctionType;
		unsigned int thisPlotID;
		ICorrection * XUnfolder;
		vector< ICorrection* > systematicUnfolders;
		vector< double > systematicOffsets, systematicWidths;
		IIndexCalculator * distributionIndices;
		TRandom3 * systematicRandom;
		string xName, priorName;
		bool finalised, normalise;
		double scaleFactor;
		vector< double > correctedDataErrors;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution, *mcRecoDistribution;
		TH2F *smearingMatrix, *covarianceMatrix;
		vector< TH1F* > systematicResults;
};

#endif
