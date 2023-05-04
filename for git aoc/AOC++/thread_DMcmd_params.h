#pragma once

#include <cstdlib>

//#include <Windows.h>
#include <cstdlib>
#include <Windows.h>

#include <MCODE.h>

using namespace MultiCODE;

typedef struct
{
	bool			bEndThread;
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;
	HANDLE			hStartServer;
	HANDLE			hStopServer;

	int				procNumber;

	DWORD			waitTimeout_ms;

	bool			bEnableDM;

	mcUINT32		loopTiming_microsec;

} ThreadParams_DMcmd;

// FUNCTION PROTOTYPES
mcINT16 StartThread_DMcmd();
mcINT16 TerminateThread_DMcmd();

int collapse_DM_buffers();
int update_DM_static_buffer();
int DM_set_to_position();

