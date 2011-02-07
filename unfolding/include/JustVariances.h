/**
  @class JustVariances

  A class to calculate the errors on the unfolded data distribution

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef JUST_VARIANCES_H
#define JUST_VARIANCES_H

#include "SmearingMatrix.h"
#include "UnfoldingMatrix.h"
#include "Distribution.h"
#include "Indices.h"
#include <vector>

using namespace std;

class JustVariances
{
	public:
		JustVariances();
		JustVariances( UnfoldingMatrix*, SmearingMatrix*, Distribution*, double, Indices* );
		~JustVariances();

		double GetVariance(int);
		double GetStandardDeviation(int);
		vector<double> GetVariances();
		vector<double> GetStandardDeviations();

	private:
		void Calculate();
		double UnfoldingCovariance( int, int, int );
		double UnfoldingCovarianceTerm( int, int );

		vector< vector<double> > fullTerms;
		vector<double> variances;
		SmearingMatrix * inputSmearing;
		UnfoldingMatrix * inputUnfolding;
		Distribution * inputDistribution;
		double inputIntegral;
		int binNumber;
};

#endif
