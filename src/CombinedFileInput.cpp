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

const double JET_PT_CUTS[] = { 30000.0, 50000.0, 100000.0, 170000.0, 320000.0, 600000.0, 1180000.0 };

CombinedFileInput::CombinedFileInput()
{
}

CombinedFileInput::CombinedFileInput( vector< string > FilePaths, vector< double > FileWeights, string InternalPath, string InputType, string Description,
		unsigned int DescriptionIndex, ObservableList * RelevanceChecker )
{
	m_sourceDescription = Description;
	m_sourceIndex = DescriptionIndex;
	m_fileWeights = FileWeights;
	m_filePaths = FilePaths;
	m_internalPath = InternalPath;
	m_inputType = InputType;
	m_currentFile = 0;
	m_rowInCurrentFile = 0;
	m_currentInput = 0;
	m_pTcuts = vector<double>( JET_PT_CUTS, JET_PT_CUTS + sizeof( JET_PT_CUTS ) / sizeof( double ) );
	m_relevanceChecker = RelevanceChecker;

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
	ChangeInputFile( 0 );
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
			ChangeInputFile( FileIndex );
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
			ChangeInputFile( FileIndex );
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
vector< double > * CombinedFileInput::GetVector( string VectorName )
{
	return m_currentInput->GetVector( VectorName );
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

void CombinedFileInput::ChangeInputFile( unsigned int NewFileIndex )
{
	m_currentFile = NewFileIndex;
	if ( m_currentInput )
	{
        	delete m_currentInput;
	}

	//Choose the type of file to open
	if ( m_inputType == NTUPLE_TYPE_STRING )
	{
		m_currentInput = new InputNtuple( m_filePaths[ m_currentFile ], m_internalPath, m_sourceDescription, m_sourceIndex, m_relevanceChecker );
	}
	else if ( m_inputType == UE_TREE_TYPE_STRING )
	{
		m_currentInput = new InputUETree( m_filePaths[ m_currentFile ], m_internalPath, m_sourceDescription, m_sourceIndex, m_relevanceChecker );
	}
	else if ( m_inputType == TRIGGER_CHOOSING_TYPE_STRING )
	{
		//if ( m_pTcuts.size() == m_filePaths.size() )
		if ( m_fileWeights[ m_currentFile ] != 1.0 )
		{
			m_currentInput = new TriggerChoosingInput( m_filePaths[ m_currentFile ], m_internalPath, m_sourceDescription, m_sourceIndex, m_relevanceChecker, m_pTcuts[ m_currentFile ] );
		}
		else
		{
			m_currentInput = new TriggerChoosingInput( m_filePaths[ m_currentFile ], m_internalPath, m_sourceDescription, m_sourceIndex, m_relevanceChecker );
		}
	}
	else
	{
		cerr << "Unrecognised input type: " << m_inputType << endl;
		exit(1);
	}
}
