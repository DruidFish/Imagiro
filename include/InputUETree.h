/**
  @class InputUETree

  Accesses a slightly less general file format I use for the underlying event analysis

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */


#ifndef INPUT_UE_TREE_H
#define INPUT_UE_TREE_H

#include "IFileInput.h"
#include <map>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "MonteCarloInformation.h"

using namespace std;

class InputUETree : public IFileInput
{
	public:
		InputUETree();
		InputUETree( string FilePath, string NtuplePath, string Description, int InputIndex );
		~InputUETree();

		//Change the Ntuple row being examined
		virtual bool ReadRow( long RowIndex );
		virtual bool ReadEvent( UInt_t EventNumber );
		virtual bool ReadEvent( UInt_t EventNumber, int FileIndex );

		//Get the standard event number and weight information
		virtual UInt_t EventNumber();
		virtual double EventWeight();

		//Get any other column value by name
		virtual double GetValue( string VariableName );

		//Get the number of rows
		virtual long NumberOfRows();
		virtual long CurrentRow();
		virtual long CurrentFile();

		//Get the description of the source
		virtual string * Description();
		virtual int DescriptionIndex();

	private:
		//Caching
		long currentRowNumber, numberOfRows;
		UInt_t currentEventNumber;
		double currentEventWeight;
		TBranch *eventNumberBranch, *eventWeightBranch;
		vector< TBranch* > branches;
		vector< double > currentValues;

		//Mapping
		map< UInt_t, long > eventNumberToRow;
		map< UInt_t, long >::iterator eventIterator;
		map< string, int > columnNameToIndex;
		map< string, int >::iterator columnIterator;

		//IO
		TFile * inputFile;
		TTree * wrappedNtuple;
		string sourceDescription;
		int sourceDescriptionIndex;
};

#endif
