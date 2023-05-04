#include <Windows.h>

#include <iostream>
#include <vector>
#include <map>
//#include "DLLfromPython.h"

using namespace std;

#define THREAD_TERMINATION_TIMEOUT	5000
#define	THREAD_CREATION_TIMEOUT		5000

map<string, HANDLE> mapEvent;

extern "C"
int register_named_event(char *chEventName, int bOpenOnly)
{
	string evt_name = string(chEventName);

	map<string, HANDLE>::iterator it = mapEvent.begin();

	for(size_t k=0;k<mapEvent.size();++k)
	{
		if ((*it).first.compare(chEventName) == 0)
		{
			cout << "Event named '" << evt_name << "' already registered. Nothing done.\n";
			return -1;
		}
	}

	HANDLE myEvent;

	if (bOpenOnly == 1)
	{
		myEvent = OpenEvent(EVENT_ALL_ACCESS,FALSE,evt_name.data());
		if (myEvent == NULL)
		{
			cout << "ERROR: failed to open event named '" << evt_name << "'.\n";
			return -1;
		}
		cout << "Event named '" << evt_name << "' has been successfully opened and registered.\n";
	}
	else
	{
		myEvent = CreateEvent(NULL,FALSE,FALSE,evt_name.data());
		if (myEvent == NULL)
		{
			cout << "ERROR: failed to create event named '" << evt_name << "'.\n";
			return -1;
		}
		cout << "Event named '" << evt_name << "' has been successfully created and registered.\n";
	}

	mapEvent.insert(std::pair<string,HANDLE>(evt_name, myEvent));

	
	return 0;
}

extern "C"
int unregister_named_event(char *chEventName)
{
	string evt_name = string(chEventName);
	
	size_t index = 0;

	bool bFound = false;

	map<string, HANDLE>::iterator it = mapEvent.begin();

	for(;it!=mapEvent.end();++it)
	{
		if ((*it).first.compare(evt_name) == 0)
		{
			bFound = true;
		}
	}

	if (bFound == false)
	{
		cout << "Event named '" << evt_name << "' not found. Could not unregister.\n";
		return -1;
	}

	CloseHandle(mapEvent[evt_name]);
	mapEvent.erase(evt_name);

	cout << "Event named '" << evt_name << "' successfully unregistered.\n";
	
	return 0;
}


extern "C"
int unregister_all()
{
	bool bError = false;

	map<string, HANDLE>::iterator it = mapEvent.begin();

	for(;it!=mapEvent.end();++it)
	{
		if (CloseHandle((*it).second) == 0)
		{
			cout << "ERROR: event named '" << (*it).first << "' could not be closed.\n";
			bError = true;
		}
		else
		{
			cout << "SUCCESS: event named '" << (*it).first << "' has been successfully closed.\n";
		}
	}

	mapEvent.clear();

	if (bError == true)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

extern "C"
int trigger_named_event(char *chEventName)
{
	HANDLE evt = INVALID_HANDLE_VALUE;

	try
	{
		evt = mapEvent.at(string(chEventName));
	}
	catch(const std::out_of_range& oor)
	{
		cout << "Error: key '" << oor.what() << "' does not exist.\n";
		return -1;
	}

	SetEvent(evt);

	cout << "Event named '" << chEventName << "' has been triggered.\n";
	
	return 0;
}


extern "C"
int chk_wait_status_named_event(char *chEventName)
{
	HANDLE evt = OpenEvent(EVENT_ALL_ACCESS,FALSE,string(chEventName).data());

	if (evt == NULL)
	{
		return -1;
	}

	int waitResult = WaitForSingleObject(evt, 0);

	CloseHandle(evt);
	
	return waitResult;
}

extern "C"
int wait_for_single_event(char *chEventName, DWORD timeout)
{
	HANDLE evt = OpenEvent(EVENT_ALL_ACCESS,FALSE,string(chEventName).data());

	if (evt == NULL)
	{
		CloseHandle(evt);
		cout << "Error: event pointer is NULL!\n";
		return -1;
	}

	int waitResult = WaitForSingleObject(evt, timeout);
	
	CloseHandle(evt);

	return waitResult;
}

extern "C"
int wait_for_multiple_events(int nEvents, bool bWaitForAllEvents, const char **chEventNames, DWORD timeout)
{
	if (nEvents <= 0)
	{
		cout << "Error: invalid 'nEvents' value (" << nEvents << ").\n";
		return -1;
	}

	vector<HANDLE> vEvent(nEvents);

	for (int k = 0; k < nEvents; ++k)
	{
		//cout << string(chEventNames[k]) << endl;
		vEvent[k] = OpenEvent(EVENT_ALL_ACCESS,FALSE,string(chEventNames[k]).data());

		if (vEvent[k] == NULL)
		{
			for (int i = 0; i <= k; ++i)
			{
				CloseHandle(vEvent[i]);
			}
			cout << "Error: event pointer indexed " << k << " (\"" << string(chEventNames[k]) << "\")is NULL!\n";
			return -1;
		}
	}

	int waitResult = WaitForMultipleObjects(DWORD(nEvents),&vEvent[0],bWaitForAllEvents,timeout);
	
	for (int k = 0; k < nEvents; ++k)
	{
		CloseHandle(vEvent[k]);
	}

	return waitResult;
}
