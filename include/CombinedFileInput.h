/**
  @class CombinedFileInput

  Combines input files on the fly without affecting their structure on disk

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */


#ifndef COMBINED_FILE_INPUT_H
#define COMBINED_FILE_INPUT_H

#include <vector>
#include <string>
#include "IFileInput.h"

using namespace std;

class CombinedFileInput : public IFileInput
{
	public:
		CombinedFileInput();
		CombinedFileInput( vector< string > FilePaths, vector< double > FileWeights, string InternalPath, string InputType, string Description, unsigned int DescriptionIndex );
		virtual ~CombinedFileInput();

		//Access a particular event, return false if the event is not found
		virtual bool ReadRow( unsigned long RowIndex, unsigned int FileIndex );
		virtual bool ReadEvent( UInt_t EventNumber, unsigned int FileIndex );

		//Get the standard event number and weight information
		virtual UInt_t EventNumber();
		virtual double EventWeight();

		//Get any other column value by name
		virtual double GetValue( string VariableName );

		//Get the number of rows
		virtual unsigned long NumberOfRows();
		virtual unsigned long CurrentRow();
		virtual unsigned int NumberOfFiles();
		virtual unsigned int CurrentFile();

		//Get the description of the source
		virtual string * Description();
		virtual unsigned int DescriptionIndex();

	private:
		//Used to instantiate files when they are needed
		IFileInput * InstantiateSingleInput( string FilePath, string InternalPath, string Type );

		string m_sourceDescription, m_inputType, m_internalPath;
		unsigned int m_currentFile, m_sourceIndex;
		unsigned long m_rowInCurrentFile;
		vector< string > m_filePaths;
		vector< double > m_fileWeights;
		vector< long > m_eventsPerFile;
		IFileInput * m_currentInput;
};

#endif
