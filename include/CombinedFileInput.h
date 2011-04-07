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
		CombinedFileInput( vector< IFileInput* > FileInputs, vector< double > FileWeights, string Description, int DescriptionIndex );
		~CombinedFileInput();

		//Access a particular event, return false if the event is not found
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
		UInt_t currentEventNumber;
		string sourceDescription;
		int currentFile, sourceIndex;
		long totalRows;
		vector< IFileInput* > fileInputs;
		vector< double > fileWeights;
};

#endif
