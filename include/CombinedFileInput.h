/**
  @class CombinedFileInput

  Combines input files on the fly without affecting their structure on disk

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */


#ifndef COMBINED_FILE_INPUT_H
#define COMBINED_FILE_INPUT_H

#include <vector>
#include "IFileInput.h"

using namespace std;

class CombinedFileInput : public IFileInput
{
	public:
		CombinedFileInput();
		CombinedFileInput( vector< IFileInput* > FileInputs, vector< double > FileWeights, string Description, unsigned int DescriptionIndex );
		~CombinedFileInput();

		//Access a particular event, return false if the event is not found
		virtual bool ReadRow( unsigned long RowIndex );
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
		UInt_t currentEventNumber;
		string sourceDescription;
		unsigned int currentFile, sourceIndex;
		unsigned long totalRows;
		vector< IFileInput* > fileInputs;
		vector< double > fileWeights;
};

#endif
