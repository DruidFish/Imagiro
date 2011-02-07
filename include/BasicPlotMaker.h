#ifndef BASIC_PLOTMAKER_H
#define BASIC_PLOTMAKER_H

#include "IPlotMaker.h"
#include "IterativeUnfolding.h"
#include "Indices.h"
#include <string>

using namespace std;

class BasicPlotMaker : public IPlotMaker
{
	public:
		BasicPlotMaker();
		BasicPlotMaker( string XVariableName, string YVariableName, string PriorName, int XBinNumber, double XMinimum, double XMaximum,
				double ScaleFactor = 1.0, int UniqueID = 0 );
		~BasicPlotMaker();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( InputNtuple * TruthInput, InputNtuple * ReconstructedInput );
		virtual void StoreMiss( InputNtuple * TruthInput );
		virtual void StoreFake( InputNtuple * ReconstructedInput );
		virtual void StoreData( InputNtuple * DataInput );

		//Do the unfolding
		virtual void Unfold( int MostIterations = 10, double ChiSquaredThreshold = 1.0, double KolmogorovThreshold = 1.0, bool WithSmoothing = false );

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
		int uniqueID;		
		IterativeUnfolding *XvsYUnfolder, *XUnfolder;
		Indices *DistributionIndices;
		string xName, yName, priorName;
		bool finalised;
		double scaleFactor;
		vector<double> correctedDataErrors;
		TH1F *correctedDistribution, *uncorrectedDistribution, *mcTruthDistribution;
		TH2F *smearingMatrix;
};

#endif
