/**
  @class SmearingMatrix

  The crucial part of the unfolding process - the matrix describing detector effects on data
  Sparse matrix version thereof

  @author Benjamin M Wynne bwynne@cern.ch
  @date 08-04-2011
 */

#ifndef SMEARING_MATRIX_H
#define SMEARING_MATRIX_H

#include <vector>
#include <map>
#include "IIndexCalculator.h"
#include "SparseMatrix.h"

using namespace std;

class SmearingMatrix : public SparseMatrix
{
	public:
		SmearingMatrix();
		SmearingMatrix( IIndexCalculator * InputIndices );
		~SmearingMatrix();

		void StoreTruthRecoPair( vector< double > Truth, vector< double > Reco, double TruthWeight = 1.0, double RecoWeight = 1.0 );
		void StoreUnreconstructedTruth( vector< double > Truth, double Weight = 1.0 );
		void StoreReconstructedFake( vector< double > Reco, double Weight = 1.0 );

		void Finalise();

                double GetEfficiency( unsigned int CauseIndex );
		double GetTruthTotal( unsigned int CauseIndex );
		double GetTotalPaired();
		double GetTotalMissed();
		double GetTotalFake();

	private:
		double totalPaired, totalMissed, totalFake;
		vector< double > normalisation, efficiencies;
		IIndexCalculator * indexCalculator;
		bool isFinalised;
};

#endif
