#ifndef X_VS_Y_NORMALISED_FOLDING_H
#define X_VS_Y_NORMALISED_FOLDING_H

#include "IPlotMaker.h"
#include "Folding.h"
#include "DataIndices.h"
#include <string>

using namespace std;

class XvsYNormalisedFolding : public IPlotMaker
{
	public:
		XvsYNormalisedFolding();
		XvsYNormalisedFolding( string XVariableName, string YVariableName, string PriorName, int XBinNumber, double XMinimum, double XMaximum,
				int YBinNumber, double YMinimum, double YMaximum, double ScaleFactor = 1.0, int UniqueID = 0 );
		~XvsYNormalisedFolding();

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
		//WARNING: this method deletes the argument object
		TH1F * Delinearise( TH1F * LinearisedDistribution );
		vector<double> DelineariseErrors( vector<double> InputSumWeightSquares );

		int uniqueID;		
		Folding *XvsYFolder, *XFolder;
		DataIndices *DistributionIndices;
		string xName, yName, priorName;
		bool finalised;
		double scaleFactor;
		vector<double> correctedInputErrors;
		TH1F *foldedDistribution, *inputDistribution, *reconstructedDistribution;
		TH2F *smearingMatrix;
};

#endif