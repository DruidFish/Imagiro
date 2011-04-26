/**
  @class UnfoldingMatrix

  The Bayesian "inverse" of the smearing matrix, used to unfold the data
  Sparse version thereof

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef UNFOLDING_MATRIX_H
#define UNFOLDING_MATRIX_H

#include "SparseMatrix.h"
#include "SmearingMatrix.h"
#include "Distribution.h"
#include "Indices.h"
#include <vector>

using namespace std;

class UnfoldingMatrix : public SparseMatrix
{
	public:
		UnfoldingMatrix();
		UnfoldingMatrix( SmearingMatrix * InputSmearing, Distribution * InputDistribution );
		~UnfoldingMatrix();
};

#endif
