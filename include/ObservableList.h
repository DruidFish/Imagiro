#ifndef OBSERVABLE_LIST_H
#define OBSERVABLE_LIST_H

//#include "MonteCarloSummaryPlotMaker.h"
#include <vector>
#include <set>
#include <string>

using namespace std;

//Forward declare MonteCarloSummaryPlotMaker os as not to get loop dependency
class MonteCarloSummaryPlotMaker;

class ObservableList
{
	public:
		ObservableList();
		ObservableList( vector< MonteCarloSummaryPlotMaker* > AllPlots );
		~ObservableList();

		bool IsInList( string TestName );

	private:
		set< string > allNames;
};

#endif
