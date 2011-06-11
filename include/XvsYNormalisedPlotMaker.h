/**
  @class XvsYNormalisedPlotMaker

  Unfolds a 2D distribution, and divides it by the normalised 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#ifndef X_VS_Y_NORMALISED_PLOTMAKER_H
#define X_VS_Y_NORMALISED_PLOTMAKER_H

#include "IUnfolder.h"
#include "StatisticsSummary.h"
#include "IterativeUnfolding.h"
#include "DataIndices.h"
#include <string>

using namespace std;

class XvsYNormalisedPlotMaker : public IUnfolder
{
	public:
		XvsYNormalisedPlotMaker();
		XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName,
				unsigned int XBinNumber, double XMinimum, double XMaximum,
				unsigned int YBinNumber, double YMinimum, double YMaximum, double ScaleFactor = 1.0 );
		~XvsYNormalisedPlotMaker();

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
		virtual vector< string > VariableNames();

	private:
		//WARNING: this method deletes the argument object
		TH1F * Delinearise( TH1F * LinearisedDistribution );
		vector< double > DelineariseErrors( vector< double > InputSumWeightSquares );

		unsigned int thisPlotID;
		IterativeUnfolding *XvsYUnfolder, *XUnfolder;
		DataIndices *DistributionIndices;
		string xName, yName, priorName;
		bool finalised;
		double scaleFactor;
		vector<double> correctedDataErrors, dagostiniErrors;
		StatisticsSummary * yValueSummary;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution, *xvsyTruthCheck, *xTruthCheck;
		TH2F *smearingMatrix, *covarianceMatrix;
		bool doPlotSmearing;
};

#endif
