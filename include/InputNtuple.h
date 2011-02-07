#ifndef INPUT_NTUPLE_H
#define INPUT_NTUPLE_H

#include <string>
#include <map>
#include <vector>
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"

using namespace std;

class InputNtuple
{
	public:
	       InputNtuple();
	       InputNtuple( string FilePath, string NtuplePath, string Description );
	       ~InputNtuple();

	       //Change the Ntuple row being examined
	       bool ReadRow( long RowIndex );
	       bool ReadEvent( float EventNumber );

	       //Get the standard event number and weight information
	       float EventNumber();
	       float EventWeight();

	       //Get any other column value by name
	       float GetValue( string VariableName );

	       //Get the number of rows
	       long NumberOfRows();
	       long CurrentRow();

	       //Get the description of the source
	       string * Description();

	private:
	       //Caching
	       long currentRowNumber, numberOfRows;
	       float currentEventNumber, currentEventWeight;
	       TBranch *eventNumberBranch, *eventWeightBranch;
	       vector< TBranch* > branches;
	       vector< float > currentValues;

	       //Mapping
	       map< float, long > eventNumberToRow;
	       map< float, long >::iterator eventIterator;
	       map< string, int > columnNameToIndex;
	       map< string, int >::iterator columnIterator;

	       //IO
	       TFile *inputFile;
	       TNtuple *wrappedNtuple;
	       string sourceDescription;
};

#endif
