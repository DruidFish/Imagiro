/**
  @class UnfoldingMatrix

  The Bayesian "inverse" of the smearing matrix, used to unfold the data

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef UNFOLDING_MATRIX_H
#define UNFOLDING_MATRIX_H

#include "SmearingMatrix.h"
#include "Distribution.h"
#include "Indices.h"
#include <vector>

using namespace std;

class UnfoldingMatrix
{
	public:
		UnfoldingMatrix();
		UnfoldingMatrix( SmearingMatrix*, Distribution*, Indices* );
		~UnfoldingMatrix();

		double GetElement( int, int );

	private:
		vector< vector<double> > matrix, covariance;
};

#endif
