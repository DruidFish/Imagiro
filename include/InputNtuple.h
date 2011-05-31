/**
  @class InputNtuple

  Loads data from a file containing a Root Ntuple, and uses hashtables to provide efficient access to that data

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-01-2011
 */


#ifndef INPUT_NTUPLE_H
#define INPUT_NTUPLE_H

#include <map>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "IFileInput.h"

using namespace std;

class InputNtuple : public IFileInput
{
	public:
		InputNtuple();
		InputNtuple( string FilePath, string NtuplePath, string Description, unsigned int InputIndex );
		~InputNtuple();

		//Change the Ntuple row being examined
		virtual bool ReadRow( unsigned long RowIndex );
		virtual bool ReadNextRow();
		virtual bool ReadEvent( UInt_t EventNumber );
		virtual bool ReadEvent( UInt_t EventNumber, unsigned int FileIndex );

		//Get the standard event number and weight information
		virtual UInt_t EventNumber();
		virtual double EventWeight();

		//Get any other column value by name
		virtual double GetValue( string VariableName );

		//Get the number of rows
		virtual unsigned long NumberOfRows();
		virtual unsigned long CurrentRow();
		virtual unsigned int CurrentFile();

		//Get the description of the source
		virtual string * Description();
		virtual unsigned int DescriptionIndex();

	private:
		//Caching
		unsigned long currentRowNumber, numberOfRows;
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
		TFile * inputFile;
		TTree * wrappedNtuple;
		string sourceDescription;
		unsigned int sourceDescriptionIndex;
};

#endif
