/**
  @class XFolding

  Applies the smearing matrix to a 1D distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */

#ifndef X_FOLDING_H
#define X_FOLDING_H

#include "IFolder.h"
#include "Folding.h"
#include "IIndexCalculator.h"
#include <string>

using namespace std;

class XFolding : public IFolder
{
	public:
		XFolding();
		XFolding( string XVariableName, string PriorName, unsigned int XBinNumber, double XMinimum, double XMaximum, double ScaleFactor = 1.0, bool Normalise = false );
		~XFolding();

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
		virtual vector<string> VariableNames();

	private:
		//To be used with Clone
		XFolding( string XVariableName, string PriorName, IIndexCalculator * DistributionIndices, double ScaleFactor = 1.0, bool Normalise = false );

		unsigned int thisPlotID;
		Folding * XFolder;
		IIndexCalculator * distributionIndices;
		string xName, priorName;
		bool finalised, normalise;
		double scaleFactor;
		vector< double > correctedInputErrors;
		TH1F *foldedDistribution, *inputDistribution, *reconstructedDistribution;
		TH2F *smearingMatrix;
};

#endif
