/**
  @interface IFileInput

  A general interface for any class that loads data from Root files

  @author Benjamin M Wynne bwynne@cern.ch
  @date 06-04-2011
 */


#ifndef I_FILE_INPUT_H
#define I_FILE_INPUT_H

#include <string>
#include "Rtypes.h"

using namespace std;

class IFileInput
{
	public:
		//Destructor
		virtual ~IFileInput()
		{
		}

		//Access a particular event, return false if the event is not found
		virtual bool ReadRow( unsigned long RowIndex, unsigned int FileIndex ) = 0;
		virtual bool ReadEvent( UInt_t EventNumber, unsigned int FileIndex ) = 0;

		//Get the standard event number and weight information
		virtual UInt_t EventNumber() = 0;
		virtual double EventWeight() = 0;

		//Get any other column value by name
		virtual double GetValue( string VariableName ) = 0;

		//Get the number of rows and files
		virtual unsigned long NumberOfRows() = 0;
		virtual unsigned long CurrentRow() = 0;
		virtual unsigned int NumberOfFiles() = 0;
		virtual unsigned int CurrentFile() = 0;

		//Get the description of the source
		virtual string * Description() = 0;
		virtual unsigned int DescriptionIndex() = 0;
};

#endif
