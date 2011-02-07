/**
  @class SmearingMatrix

  The crucial part of the unfolding process - the matrix describing detector effects on data

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef SMEARING_MATRIX_H
#define SMEARING_MATRIX_H

#include <vector>
#include <string>
#include "Indices.h"
#include "TH2F.h"

using namespace std;

class SmearingMatrix
{
	public:
		SmearingMatrix();
		SmearingMatrix( Indices* );
		~SmearingMatrix();

		void StoreTruthRecoPair( vector<double> Truth, vector<double> Reco, double Weight = 1.0 );
		void StoreUnreconstructedTruth( vector<double> Truth, double Weight = 1.0 );
		void StoreReconstructedFake( vector<double> Reco, double Weight = 1.0 );
		void StoreUnnormalisedMatrix( TH2F* );

		void Finalise();

		double GetElement( int, int );
		double GetEfficiency( int );
		double GetCovarianceElement( int, int, int );
		TH2F * MakeRootHistogram( string, string, bool MakeNormalised = true );

	private:
		vector< vector<double> > matrix;
		//vector< vector< vector<double> > > smearingCovariance;
		vector<double> efficiencies, normalisation;
		Indices * indexCalculator;
		bool isFinalised;
};

#endif
