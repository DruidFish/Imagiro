#ifndef COMPARISON_H
#define COMPARISON_H

#include "Distribution.h"
#include <string>

using namespace std;

class Comparison
{
	public:
		Comparison();
		Comparison( string Name, int UniqueID );
		~Comparison();

		void CompareDistributions( Distribution * FirstInput, Distribution * SecondInput, double & ChiSquared, double & Kolmogorov, bool Normalised, bool IsClosureTest = false );

	private:
		string name;
		int internalID, uniqueID;

};

#endif
