/**
  @class CovarianceMatrix

  A class for calculating the covariance matrix for the unfolded data distribution
  Depreciated now since it's really slow, but it might work

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef COVARIANCE_MATRIX_H
#define COVARIANCE_MATRIX_H

#include "SmearingMatrix.h"
#include "UnfoldingMatrix.h"
#include "Distribution.h"
#include "Indices.h"
#include <vector>

using namespace std;

class CovarianceMatrix
{
	public:
		CovarianceMatrix();
		CovarianceMatrix( UnfoldingMatrix*, SmearingMatrix*, Distribution*, double, Indices* );
		~CovarianceMatrix();

		double GetElement( int, int );

	private:
		double UnfoldingCovariance( int, int, int, int );
		double UnfoldingCovarianceTerm( int, int, int, int );
		void MakeUnfoldingCovarianceTerm();

		vector< vector< vector< vector<double> > > > unfoldingCovarianceTerm;
		vector< vector<double> > matrix;
		SmearingMatrix * inputSmearing;
		UnfoldingMatrix * inputUnfolding;
		int binNumber;
};

#endif
