/**
  @class XvsYNormalisedPlotMaker

  Unfolds a 2D distribution, and divides it by the normalised 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#ifndef X_VS_Y_NORMALISED_PLOTMAKER_H
#define X_VS_Y_NORMALISED_PLOTMAKER_H

#include "IPlotMaker.h"
#include "StatisticsSummary.h"
#include "ICorrection.h"
#include "IIndexCalculator.h"
#include "TProfile.h"

using namespace std;

class XvsYNormalisedPlotMaker : public IPlotMaker
{
	public:
		XvsYNormalisedPlotMaker();
		XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
				unsigned int XBinNumber, double XMinimum, double XMaximum,
				unsigned int YBinNumber, double YMinimum, double YMaximum, int CorrectionMode = 2, double ScaleFactor = 1.0 );
		XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
				vector< double > XBinLowEdges, vector< double > YBinLowEdges, int CorrectionMode = 2, double ScaleFactor = 1.0 );

		virtual ~XvsYNormalisedPlotMaker();

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
		virtual XvsYNormalisedPlotMaker * Clone( string NewPriorName );

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
		//To be used only with Clone
		XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
				IIndexCalculator * DistributionIndices, int CorrectionMode, unsigned int OriginalID, double ScaleFactor = 1.0 );

		//Instantiate an object to correct the data
		ICorrection * MakeCorrector( int CorrectionMode, IIndexCalculator * CorrectionIndices, string CorrectionName, unsigned int CorrectionID );

		//WARNING: this method deletes the argument object
		TH1F * MakeProfile( TH1F * LinearisedDistribution );
		vector< double > DelineariseErrors( vector< double > InputSumWeightSquares );

		TProfile * simpleDataProfile; 
		int correctionType;
		unsigned int thisPlotID, xBinNumber, yBinNumber;
		ICorrection *XvsYUnfolder;
		IIndexCalculator *distributionIndices;
		string xName, yName, priorName;
		bool finalised, doPlotSmearing;
		double scaleFactor, xMinimum, yMinimum, xMaximum, yMaximum;
		vector< double > correctedDataErrors;
		StatisticsSummary * yValueSummary;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution, *xvsyTruthCheck, *xTruthCheck;
		TH2F *smearingMatrix, *covarianceMatrix;
};

#endif
