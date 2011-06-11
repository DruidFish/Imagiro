/**
  @class XvsYNormalisedFolding

  Folds a 2D distribution, and divides it by the folded 1D distribution of the X-axis variable (giving value-of-Y-per-event).
  Includes error checking on the discretisation of the Y variable

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#ifndef X_VS_Y_NORMALISED_FOLDING_H
#define X_VS_Y_NORMALISED_FOLDING_H

#include "IFolder.h"
#include "StatisticsSummary.h"
#include "Folding.h"
#include "DataIndices.h"
#include <string>

using namespace std;

class XvsYNormalisedFolding : public IFolder
{
	public:
		XvsYNormalisedFolding();
		XvsYNormalisedFolding( string XVariableName, string YVariableName, string PriorName,
				unsigned int XBinNumber, double XMinimum, double XMaximum,
				unsigned int YBinNumber, double YMinimum, double YMaximum, double ScaleFactor = 1.0 );
		~XvsYNormalisedFolding();

		//Take input values from ntuples
		//To reduce file access, the appropriate row must already be in memory, the method does not change row
		virtual void StoreMatch( IFileInput * TruthInput, IFileInput * ReconstructedInput );
		virtual void StoreMiss( IFileInput * TruthInput );
		virtual void StoreFake( IFileInput * ReconstructedInput );
		virtual void StoreData( IFileInput * DataInput );

		//Do the folding
		virtual void Fold();

		//Do a closure test
                virtual bool ClosureTest();

		//Return some plots
		virtual TH1F * CorrectedHistogram();
		virtual TH1F * UncorrectedHistogram();
		virtual TH1F * MCTruthHistogram();
		virtual TH2F * SmearingMatrix();

		//Copy the object
		virtual IFolder * Clone( string NewPriorName );

		//General info
		virtual string Description( bool WithSpaces );
		virtual string PriorName();

		//Error info for corrected distribution
		virtual vector< double > CorrectedErrors();

		//Return the names of the variables involved
		virtual vector< string > VariableNames();

	private:
		//WARNING: this method deletes the argument object
		TH1F * Delinearise( TH1F * LinearisedDistribution );
		vector< double > DelineariseErrors( vector< double > InputSumWeightSquares );

		unsigned int thisPlotID;
		Folding *XvsYFolder, *XFolder;
		DataIndices *DistributionIndices;
		string xName, yName, priorName;
		bool finalised;
		double scaleFactor;
		vector< double > correctedInputErrors;
		StatisticsSummary * yValueSummary;
		TH1F *foldedDistribution, *inputDistribution, *reconstructedDistribution, *xvsyReconstructionCheck, *xReconstructionCheck;
		TH2F *smearingMatrix;
};

#endif
