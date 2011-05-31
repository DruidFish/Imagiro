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

//CombinedFileInput::CombinedFileInput( vector< IFileInput* > FileInputs, vector< double > FileWeights, string Description, unsigned int DescriptionIndex )
CombinedFileInput::CombinedFileInput( vector< string > FilePaths, vector< double > FileWeights, string InternalPath, string InputType, string Description, unsigned int DescriptionIndex )
{
	m_sourceDescription = Description;
	m_sourceIndex = DescriptionIndex;
	m_fileWeights = FileWeights;
	m_filePaths = FilePaths;
	m_internalPath = InternalPath;
	m_inputType = InputType;
	m_totalRows = 0;
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

	cout << "Initialisation loading file 0 " << FilePaths[ 0 ] << endl;
	m_currentInput = InstantiateSingleInput( FilePaths[ 0 ], InternalPath, InputType );
	m_rowInCurrentFile = m_currentInput->CurrentRow();
	m_eventsPerFile = vector< long >( FilePaths.size(), 0 );
	m_eventsPerFile[ 0 ] = m_currentInput->NumberOfRows();

	/*else
	  {
	//Set up each slice input
	fileInputs = FileInputs;
	for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
	{
	currentFile = fileIndex;
	totalRows += fileInputs[ fileIndex ]->NumberOfRows();
	}
	}*/
}

CombinedFileInput::~CombinedFileInput()
{
	delete m_currentInput;

	/*for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
	  {
	  delete fileInputs[ fileIndex ];
	  }*/
}

//Access a particular event, return false if the event is not found
bool CombinedFileInput::ReadRow( unsigned long RowNumber )
{
	if ( RowNumber == 0 )
	{
		if ( m_currentFile == 0 )
		{
			m_rowInCurrentFile = 0;
		}
		else
		{
			m_currentFile = 0;
			m_rowInCurrentFile = 0;
			delete m_currentInput;
			cout << "ReadRow loading file " << m_currentFile << " " << m_filePaths[ m_currentFile ] << endl;
			m_currentInput = InstantiateSingleInput( m_filePaths[ m_currentFile ], m_internalPath, m_inputType );
			m_eventsPerFile[ m_currentFile ] = m_currentInput->NumberOfRows();
		}

		return m_currentInput->ReadRow( RowNumber );
	}
	else
	{
		cerr << "Trying to use CombinedFileInput::ReadRow for index other than 0" << endl;
		exit(1);
	}
}

//Access the next row
bool CombinedFileInput::ReadNextRow()
{
	if ( m_rowInCurrentFile < m_currentInput->NumberOfRows() - 1 )
	{
		//Increment position in current file
		m_rowInCurrentFile++;
		return m_currentInput->ReadRow( m_rowInCurrentFile );
	}
	else
	{
		if ( m_currentFile < m_filePaths.size() - 1 )
		{
			//Open the next file
			m_rowInCurrentFile = 0;
			m_currentFile++;
			delete m_currentInput;
			cout << "ReadNextRow loading file " << m_currentFile << " " << m_filePaths[ m_currentFile ] << endl;
			m_currentInput = InstantiateSingleInput( m_filePaths[ m_currentFile ], m_internalPath, m_inputType );
			m_eventsPerFile[ m_currentFile ] = m_currentInput->NumberOfRows();
			return m_currentInput->ReadRow( m_rowInCurrentFile );
		}
		else
		{
			//No more files
			return false;
		}
	}
}

//Access a particular event, return false if the event is not found
bool CombinedFileInput::ReadEvent( UInt_t EventNumber )
{
	//Shouldn't use this any more...
	cerr << "Trying to use CombinedFileInput::ReadRow - don't!" << endl;
	exit(1);

	/*bool found = false;

	//Look for the event in each slice
	for ( unsigned int fileIndex = 0; fileIndex < fileInputs.size(); fileIndex++ )
	{
	if ( fileInputs[ fileIndex ]->ReadEvent( EventNumber ) )
	{
	//Check there are no duplicate event numbers - should use the other ReadEvent method if so
	if ( found )
	{
	cerr << "Duplicate event numbers!" << endl;
	exit(1);
	}
	else
	{
	currentFile = fileIndex;
	found = true;
	}
	}
	}

	return found;*/
}

bool CombinedFileInput::ReadEvent( UInt_t EventNumber, unsigned int FileIndex )
{
	//See if this file is already loaded
	if ( m_currentFile != FileIndex )
	{
		if ( FileIndex >= m_filePaths.size() )
		{
			return false;
		}
		else
		{
			//Change the current input
			m_currentFile = FileIndex;
			delete m_currentInput;
			cout << "ReadEvent loading file " << m_currentFile << " " << m_filePaths[ m_currentFile ] << endl;
			m_currentInput = InstantiateSingleInput( m_filePaths[ FileIndex ], m_internalPath, m_inputType );
			m_eventsPerFile[ m_currentFile ] = m_currentInput->NumberOfRows();
		}
	}

	bool result = m_currentInput->ReadEvent( EventNumber );
	m_rowInCurrentFile = m_currentInput->CurrentRow();
	return result;

	//Constrain the search to a particular file to avoid event number conflicts
	//currentFile = FileIndex;
	//return fileInputs[ FileIndex ]->ReadEvent( EventNumber );
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
	return m_totalRows;
}

unsigned long CombinedFileInput::CurrentRow()
{
	long currentRow = 0;

	//Add up the rows in earlier slices
	for ( unsigned int fileIndex = 0; fileIndex < m_currentFile; fileIndex++ )
	{
		//cout << m_eventsPerFile[ fileIndex ] << ", ";
		currentRow += m_eventsPerFile[ fileIndex ];
	}

	//Then add the row of the active slice
	//cout << m_rowInCurrentFile << endl;
	currentRow += m_rowInCurrentFile;

	return currentRow;
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
