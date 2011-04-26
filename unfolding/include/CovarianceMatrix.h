/**
  @class CovarianceMatrix

  An attempt to implement an efficient version of D'Agostini's error propagation using sparse matrices

  @author Benjamin M Wynne bwynne@cern.ch
  @date 08-04-2011
  */

#ifndef COVARIANCE_MATRIX_H
#define COVARIANCE_MATRIX_H

#include "SmearingCovariance.h"
#include "SparseMatrix.h"
#include "UnfoldingMatrix.h"
#include "SmearingMatrix.h"
#include "Distribution.h"
#include <vector>
#include <ctime>

using namespace std;

class CovarianceMatrix : public SparseMatrix
{
	public:
		CovarianceMatrix();
		CovarianceMatrix( UnfoldingMatrix * InputUnfolding, SmearingMatrix * InputSmearing, Distribution * DataDistribution, double CorrectedSum, bool JustVariance = false );
		~CovarianceMatrix();

	private:
		void CovarianceCalculation( int I, int J, int K, int L, double unfoldingProductTimesDataI, double dataI, double dataJ );

		time_t timeNow;
		double correctedSum;
		SmearingMatrix * inputSmearing;
		UnfoldingMatrix * inputUnfolding;
		SmearingCovariance * rsuMatrix;
};

#endif
