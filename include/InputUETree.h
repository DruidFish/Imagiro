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
		InputUETree( string FilePath, string NtuplePath, string Description, unsigned int InputIndex,
				string CutName = "NoCut", double CutMinimum = 0.0, double CutMaximum = 0.0 );
		virtual ~InputUETree();

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
		unsigned long currentRowNumber;
		UInt_t currentEventNumber;
		double currentEventWeight;
		TBranch *eventNumberBranch, *eventWeightBranch;
		vector< TBranch* > branches;
		vector< double > currentValues;

		//Mapping
		map< UInt_t, unsigned long > eventNumberToRow;
		map< UInt_t, unsigned long >::iterator eventIterator;
		map< string, unsigned int > columnNameToIndex;
		map< string, unsigned int >::iterator columnIterator;
		vector< UInt_t > validRows;

		//IO
		TFile * inputFile;
		TTree * wrappedNtuple;
		string sourceDescription;
		unsigned int sourceDescriptionIndex;
};

#endif
