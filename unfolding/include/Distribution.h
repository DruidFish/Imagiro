/**
  @class Distribution

  A 1D data histogram / probability distribution. Contains the folding and unfolding methods

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
		Distribution( Indices * InputIndices );
		Distribution( vector< TH1F* > InputDistributions, Indices * InputIndices );
		Distribution( Distribution * DataDistribution, SmearingMatrix * Smearing, Distribution * PriorDistribution );
		Distribution( Distribution * InputDistribution, SmearingMatrix * Smearing );
		~Distribution();

		void StoreEvent( vector<double> Value, double Weight = 1.0 );
		void StoreBadEvent( double Weight = 1.0 );
		void SetBadBin( double Ratio );

		double GetBinNumber( int BinIndex );
		double GetBinProbability( int BinIndex );

		TH1F * MakeRootHistogram( string Name, string Title, bool MakeNormalised = false, bool WithBadBin = false );

		//Try and account for statistical fluctuations with a moving-average smearing
		void Smooth( int SideBinNumber = 1 );

	protected:
		JustVariances * errorCalculator;
		Indices * indexCalculator;
		vector<double> binValues;
		double integral;
};

#endif
