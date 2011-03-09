/**
        @class StatisticsSummary

        A class that maintains a statistical summary of events passed to it

        @author Benjamin M Wynne bwynne@cern.ch
	@date 09-03-2011
*/

#ifndef STATISTICS_SUMMARY_H
#define STATISTICS_SUMMARY_H

#include <vector>

using namespace std;

class StatisticsSummary
{
	public:
		StatisticsSummary();
		~StatisticsSummary();

		void StoreEvent( double Value, double Weight );

		double Mean();
		double Variance();
		double OptimumBinWidth();
		int OptimumBinNumber( double Minimum, double Maximum );
		double Maximum();
		double Minimum();

	private:
		double currentMaximum, currentMinimum, meanNumerator, meanDenominator, meanSquaredNumerator, meanSquaredDenominator;
		bool freshStart;
};

#endif
