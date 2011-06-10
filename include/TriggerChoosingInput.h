/**
  @class TriggerChoosingInput

  UE ANALYSIS SPECIFIC HARD CODE
  An interface to pick events from different files, according to LeadJetPt

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-05-2011
 */


#ifndef TRIGGER_CHOOSING_INPUT_H
#define TRIGGER_CHOOSING_INPUT_H

#include "IFileInput.h"
#include <string>
#include <vector>
#include "Rtypes.h"

using namespace std;

class TriggerChoosingInput : public IFileInput
{
	public:
		TriggerChoosingInput();
		TriggerChoosingInput( string FilePath, string NtuplePath, string Description, unsigned int InputIndex, double JetPtMax = 1.0E20 );
		virtual ~TriggerChoosingInput();

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
		string ReplaceString( string Input, string FindString, string ReplaceString );

		//Caching
		unsigned int currentFileNumber;
		unsigned long totalRows, currentRowNumber;

		//IO
		vector< IFileInput* > triggerInputs;
		string sourceDescription;
		unsigned int sourceDescriptionIndex;
};

#endif
