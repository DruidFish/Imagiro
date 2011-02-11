#ifndef X_VS_Y_NORMALISED_PLOTMAKER_H
#define X_VS_Y_NORMALISED_PLOTMAKER_H

#include "IPlotMaker.h"
#include "IterativeUnfolding.h"
#include "DataIndices.h"
#include <string>

using namespace std;

class XvsYNormalisedPlotMaker : public IPlotMaker
{
	public:
		XvsYNormalisedPlotMaker();
		XvsYNormalisedPlotMaker( string XVariableName, string YVariableName, string PriorName, int XBinNumber, double XMinimum, double XMaximum,
				int YBinNumber, double YMinimum, double YMaximum, double ScaleFactor = 1.0, int UniqueID = 0 );
		~XvsYNormalisedPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput );
		virtual void StoreMiss( InputNtuple * TruthInput );
		virtual void StoreFake( InputNtuple * ReconstructedInput );
		virtual void StoreData( InputNtuple * DataInput );

		//Do the unfolding
		virtual void Unfold( int MostIterations = 10, double ChiSquaredThreshold = 1.0, double KolmogorovThreshold = 1.0, bool WithSmoothing = false );

		//Make a cross-check with MC
                virtual int MonteCarloCrossCheck( TH1F * ReferencePlot, double & ChiSquaredThreshold, double & KolmogorovThreshold, bool WithSmoothing = false );

		//Return some plots
		virtual TH1F * CorrectedDistribution();
		virtual TH1F * UncorrectedDistribution();
		virtual TH1F * MCTruthDistribution();
		virtual TH2F * SmearingMatrix();

		//Copy the object
		virtual IPlotMaker * Clone( string NewPriorName );

		//General info
		virtual string Description( bool WithSpaces );
		virtual string PriorName();

                //Error info for corrected distribution
                virtual vector<double> CorrectedErrors();

	private:
		//WARNING: this method deletes the argument object
		TH1F * Delinearise( TH1F * LinearisedDistribution );
		vector<double> DelineariseErrors( vector<double> InputSumWeightSquares );

		int uniqueID;		
		IterativeUnfolding *XvsYUnfolder, *XUnfolder;
		DataIndices *DistributionIndices;
		string xName, yName, priorName;
		bool finalised;
		double scaleFactor;
		vector<double> correctedDataErrors;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution, *xvsyTruthCheck, *xTruthCheck;
		TH2F *smearingMatrix;
};

#endif
