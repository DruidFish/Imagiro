/**
  @class CombinedFileInput

  Combines input files on the fly without affecting their structure on disk

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */

#include "CombinedFileInput.h"
#include "InputNtuple.h"
#include "InputUETree.h"
#include "TriggerChoosingInput.h"
#include <iostream>
#include <cstdlib>

const string NTUPLE_TYPE_STRING = "InputNtuple";
const string UE_TREE_TYPE_STRING = "InputUETree";
const string TRIGGER_CHOOSING_TYPE_STRING = "TriggerChoosingInput";

CombinedFileInput::CombinedFileInput()
{
}

CombinedFileInput::CombinedFileInput( vector< string > FilePaths, vector< double > FileWeights, string InternalPath, string InputType, string Description, unsigned int DescriptionIndex )
{
	m_sourceDescription = Description;
	m_sourceIndex = DescriptionIndex;
	m_fileWeights = FileWeights;
	m_filePaths = FilePaths;
	m_internalPath = InternalPath;
	m_inputType = InputType;
	m_currentFile = 0;
	m_rowInCurrentFile = 0;

	//Check the input
	if ( FilePaths.size() < 1 )
	{
		cerr << "CombinedFileInput requires at least one input file" << endl;
		exit(1);
	}
	else if ( FilePaths.size() != FileWeights.size() )
	{
		cerr << "Mismatch in number of input files and number of file weights" << endl;
		exit(1);
	}

	//Load the first file
	m_currentInput = InstantiateSingleInput( FilePaths[ 0 ], InternalPath, InputType );
}

CombinedFileInput::~CombinedFileInput()
{
	m_filePaths.clear();
	m_fileWeights.clear();
	delete m_currentInput;
}

//Access a particular event, return false if the event is not found
bool CombinedFileInput::ReadRow( unsigned long RowNumber, unsigned int FileIndex )
{
	//Check we're asking for a valid file
	if ( FileIndex < m_filePaths.size() )
	{
		//See if we've already got the file loaded
		if ( FileIndex != m_currentFile )
		{
			//Load the new file
			m_currentFile = FileIndex;
			delete m_currentInput;
			m_currentInput = InstantiateSingleInput( m_filePaths[ m_currentFile ], m_internalPath, m_inputType );
		}

		//Read the row in the file
		return m_currentInput->ReadRow( RowNumber, 0 );
	}
	else
	{
		cerr << "Requested invalid file " << FileIndex << " of " << m_filePaths.size() << endl;
		exit(1);
	}
}

//Access a particular event, return false if the event is not found
bool CombinedFileInput::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
{
	//Check we're asking for a valid file
	if ( FileIndex < m_filePaths.size() )
	{
		//See if we've already got the file loaded
		if ( FileIndex != m_currentFile )
		{
			//Load the new file
			m_currentFile = FileIndex;
			delete m_currentInput;
			m_currentInput = InstantiateSingleInput( m_filePaths[ m_currentFile ], m_internalPath, m_inputType );
		}

		//Search for the event in the file
		return m_currentInput->ReadEvent( EventNumber, 0 );

	}
	else
	{
		cerr << "Requested invalid file " << FileIndex << " of " << m_filePaths.size() << endl;
		exit(1);
	}
}

//Get the standard event number and weight information
UInt_t CombinedFileInput::EventNumber()
{
	return m_currentInput->EventNumber();
}

double CombinedFileInput::EventWeight()
{
	return m_currentInput->EventWeight() * m_fileWeights[ m_currentFile ];
}

//Get any other column value by name
double CombinedFileInput::GetValue( string VariableName )
{
	return m_currentInput->GetValue( VariableName );
}

//Get the number of rows
unsigned long CombinedFileInput::NumberOfRows()
{
	return m_currentInput->NumberOfRows();
}

unsigned long CombinedFileInput::CurrentRow()
{
	return m_currentInput->CurrentRow();
}
unsigned int CombinedFileInput::NumberOfFiles()
{
	return m_filePaths.size();
}
unsigned int CombinedFileInput::CurrentFile()
{
	return m_currentFile;
}

//Get the description of the source
string * CombinedFileInput::Description()
{
	return &m_sourceDescription;
}

unsigned int CombinedFileInput::DescriptionIndex()
{
	return m_sourceIndex;
}

IFileInput * CombinedFileInput::InstantiateSingleInput( string FilePath, string InternalPath, string Type )
{
	if ( Type == NTUPLE_TYPE_STRING )
	{
		return new InputNtuple( FilePath, InternalPath, m_sourceDescription, m_sourceIndex );
	}
	else if ( Type == UE_TREE_TYPE_STRING )
	{
		return new InputUETree( FilePath, InternalPath, m_sourceDescription, m_sourceIndex );
	}
	else if ( Type == TRIGGER_CHOOSING_TYPE_STRING )
	{
		return new TriggerChoosingInput( FilePath, InternalPath, m_sourceDescription, m_sourceIndex );
	}
	else
	{
		cerr << "Unrecognised input type: " << Type << endl;
		exit(1);
	}
}
