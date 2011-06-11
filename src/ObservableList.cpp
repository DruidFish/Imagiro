#include "ObservableList.h"
#include "MonteCarloSummaryPlotMaker.h"

//Default constructor
ObservableList::ObservableList()
{
}

ObservableList::ObservableList( vector< MonteCarloSummaryPlotMaker* > AllPlots )
{
	for ( unsigned int plotIndex = 0; plotIndex < AllPlots.size(); plotIndex++ )
	{
		vector< string > variableNames = AllPlots[ plotIndex ]->VariableNames();
		allNames.insert( variableNames.begin(), variableNames.end() );
	}
}

ObservableList::~ObservableList()
{
	allNames.clear();
}

bool ObservableList::IsInList( string TestName )
{
	return ( allNames.find( TestName ) != allNames.end() );
}
