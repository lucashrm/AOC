#pragma once

#include <string>
#include <map>

#include <MCODE.h>
#include <MCODE_SciLib.h>

using namespace MultiCODE;

typedef struct
{
	// THREAD CREATION/TERMINATION EVENTS
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;

	// TRIGGERING EVENTS
	HANDLE			hMissedFrame;
	HANDLE			hStatsUpdate;

	HANDLE			hIn_AOCGeneralTrigger;			// GENERAL AOC TRIGGER
	HANDLE			hIn_PupMask_available;			// TELLS THAT THE PUPIL MASK HAS BEE UPDATED
	HANDLE			hIn_ZER_cal_ready;				// TELLS THAT ZERNIKE CALIBRATION PARAMETERS HAVE BEEN UPDATED (AND READY FOR ACTUAL CALIBRATION)
	HANDLE			hIn_ZON_cal_ready;				// TELLS THAT ZONAL CALIBRATION PARAMETERS HAVE BEEN UPDATED (AND READY FOR ACTUAL CALIBRATION)
	HANDLE			hIn_DM_static_map_updated;		// TELLS THAT THE DM STATIC MAP HAS BEEN UPDATED
	HANDLE			hIn_RefCentroidsSet;			// TELLS THAT THE CURRENT CENTROID POSITIONS MUST BE SET AS THE REFERENCE
	
	// SIGNALING TO EXTERNAL PROCESSES
	HANDLE			hOut_OK;
	HANDLE			hOut_CalibDone;

	// THREAD VARIABLES
	bool			bEndThread;
	DWORD			WaitTimeout;

	// SPECIAL HANDLES
	std::map<std::string, HANDLE>	mapHandle;
} ThreadParams_EvtWatch;


mcINT16 StartThread_EvtWatch();

mcINT16 TerminateThread_EvtWatch();

static mcUINT32 Thread_EvtWatch(void *pParam);
