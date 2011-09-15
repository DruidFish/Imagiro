/**
  @class Comparison

  An extremely simple class to return the Chi2 and Kolmogorov-Smirnoff comparison values betweeen two histograms
  The histogram naming scheme prevents the tedious complaints from Root about making multiple objects with the same name

  @author Benjamin M Wynne bwynne@cern.ch
  @date 08-03-2011
 */

#ifndef COMPARISON_H
#define COMPARISON_H

#include "IIndexCalculator.h"
#include "Distribution.h"
#include "TH1F.h"
#include <string>

using namespace std;

class Comparison
{
	public:
		Comparison();
		Comparison( string Name, int UniqueID );
		~Comparison();

		void CompareDistributions( Distribution * FirstInput, Distribution * SecondInput, double & ChiSquared, double & Kolmogorov, bool IsClosureTest = false );
		void DelineariseAndCompare( Distribution * FirstInput, Distribution * SecondInput, double & ChiSquared, double & Kolmogorov, IIndexCalculator * InputIndices );

	private:
		TH1F * MakeProfile( TH1F * LinearisedDistribution, IIndexCalculator * InputIndices );

		string name;
		int internalID, uniqueID;

};

#endif
