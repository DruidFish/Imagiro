/**
  @class Distribution

  A 1D data histogram / probability distribution. Contains the unfolding method

  @author Benjamin M Wynne bwynne@cern.ch
  @date 17-06-2010
 */

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
#include <string>
#include "Indices.h"
#include "SmearingMatrix.h"
#include "TH1F.h"

using namespace std;

//Forward declaration of JustVariances, to avoid cyclic dependency
class JustVariances;

class Distribution
{
	public:
		Distribution();
		Distribution( Indices* );
		Distribution( vector< TH1F* >, Indices* );
		Distribution( Distribution*, SmearingMatrix*, Distribution*, Indices* );
		~Distribution();

		void StoreEvent( vector<double> Value, double Weight = 1.0 );

		double GetBinNumber(int);
		double GetBinProbability(int);
		TH1F * MakeRootHistogram( string, string, bool WithErrors = false, bool MakeNormalised = false );
		void Smooth( int SideBinNumber = 1 );

	protected:
		JustVariances * errorCalculator;
		Indices * indexCalculator;
		vector<double> binValues;
		double integral;
};

#endif
