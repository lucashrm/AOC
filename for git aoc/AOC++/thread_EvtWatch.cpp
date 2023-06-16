#include <iostream>
#include <cstdlib>

#include <MultiCODE.h>
#include <MCODE.h>

#include <MCODE_SystemUtilities.h>
#include <MCODE_StdFunc.h>
#include <MCODE_Export.h>
#include <MCODE_Export_FITS.h>

#include <iostream>
#include <fstream>
#include <io.h>

#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fitsio.h>

//#include <synchapi.h>

#include <ATMCD32D.H>

#undef _TIMESPEC_DEFINED
#include <pthread.h>

#include "ImageStruct.h"
#include "ImageCreate.h"

#include "thread_EvtWatch_params.h"
#include "thread_iXon_params.h"
#include "thread_calcSH_params.h"
#include "thread_DMcmd_params.h"

extern MultiCODE::MCODE_LogManager	*Log;
extern string						LogMsg;

extern ThreadParams_EvtWatch		TP_EvtWatch;
extern DWORD						Thread_EvtWatch_ID;
extern HANDLE						hThread_EvtWatch;

extern ThreadParams_iXon			TP_iXon;
extern DWORD						Thread_iXon_ID;
extern HANDLE						hThread_iXon;

extern ThreadParams_calcSH			TP_calcSH;
extern DWORD						Thread_calcSH_ID;
extern HANDLE						hThread_calcSH;

extern IMAGE				*SHM_imarray;			// SHARED MEMORY INSTANT IXON IMAGE
extern IMAGE				*SHM_Slopes;			// SHARED MEMORY N FRAMES HISTORY DATA (DISPLAY RATE)
extern IMAGE				*SH_RefCent;			// SHARED MEMORY REFERENCE CENTROIDS DATA
extern IMAGE				*SHM_DM_slope_mask;		// SHARED MEMORY PUPIL MASK DATA
extern IMAGE				*SHM_DM_state;			// SHARED MEMORY DM STATE
extern IMAGE				*SHM_ModesInfo;			// SHARED MEMORY MODE INFORMATION
extern vector<IMAGE*>		vSHM_DM_BUFFERS;		// SHARED MEMORY

extern mcINT16						ApplyRefCentroids();

using namespace std;

mcINT16 StartThread_EvtWatch()
{
	string threadName = "EvtWatch";

	// THREAD CREATION/TERMINATION EVENTS
	TP_EvtWatch.hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_EvtWatch.hTerminateThread = CreateEvent(NULL,FALSE,FALSE,NULL);

	// TRIGGERING EVENTS
	TP_EvtWatch.hStatsUpdate = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_EvtWatch.hMissedFrame = CreateEvent(NULL,FALSE,FALSE,NULL);

	TP_EvtWatch.hIn_AOCGeneralTrigger = CreateEvent(NULL,FALSE,FALSE,"_AOC_GENERAL_TRIGGER");
	TP_EvtWatch.hIn_PupMask_available = CreateEvent(NULL,FALSE,FALSE,"_AOC_PUPMASK_AVAILABLE");
	TP_EvtWatch.hIn_ZER_cal_ready = CreateEvent(NULL,FALSE,FALSE,"_AOC_ZER_CAL_READY");
	TP_EvtWatch.hIn_ZON_cal_ready = CreateEvent(NULL,FALSE,FALSE,"_AOC_ZON_CAL_READY");
	TP_EvtWatch.hIn_DM_static_map_updated = CreateEvent(NULL,FALSE,FALSE,"_AOC_DM_STATIC_MAP_UPDATED");
	TP_EvtWatch.hIn_RefCentroidsSet = CreateEvent(NULL,FALSE,FALSE,"_AOC_REF_CENTROIDS_SET");

	TP_EvtWatch.WaitTimeout = INFINITE;
	
	// SIGNALING TO EXTERNAL PROCESSES
	TP_EvtWatch.hOut_OK = CreateEvent(NULL,FALSE,FALSE,"_AOC_OK");
	TP_EvtWatch.hOut_CalibDone = CreateEvent(NULL,FALSE,FALSE,"_AOC_CALIB_DONE");


	// THREAD VARIABLE INITIALIZATION

	hThread_EvtWatch = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_EvtWatch,&TP_EvtWatch,THREAD_PRIORITY_NORMAL,&Thread_EvtWatch_ID);

	DWORD WaitResult = WaitForSingleObject(TP_EvtWatch.hThreadCreated,THREAD_CREATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		LogMsg = string_format("[%s] thread %s created successuflly.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SUCCESS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		break;
	case WAIT_TIMEOUT:
		LogMsg = string_format("[%s] %s thread creation TIMEOUT.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_ABANDONED:
		LogMsg = string_format("[%s] %s thread creation ABANDONED.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_FAILED:
		LogMsg = string_format("[%s] %s thread creation FAILED.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	}

	return SUCCESS;
}

mcINT16 TerminateThread_EvtWatch()
{
	string threadName = "EvtWatch";

	SetEvent(TP_EvtWatch.hTerminateThread);

	DWORD WaitResult = WaitForSingleObject(hThread_EvtWatch,THREAD_TERMINATION_TIMEOUT);

	switch(WaitResult)
	{
	case WAIT_OBJECT_0:
		LogMsg = string_format("[%s] thread %s terminated successuflly.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SUCCESS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		break;
	case WAIT_TIMEOUT:
		LogMsg = string_format("[%s] %s thread termination TIMEOUT.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_ABANDONED:
		LogMsg = string_format("[%s] %s thread termination ABANDONED.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_FAILED:
		LogMsg = string_format("[%s] %s thread termination FAILED.",__FUNCTION__,threadName.data());
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	}

	return SUCCESS;
}

mcUINT32 Thread_EvtWatch(void *pParam)
{
	ThreadParams_EvtWatch *param = (ThreadParams_EvtWatch*)pParam;

	MCODE_Export_FITS FITS;

	DWORD WaitResult;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hMissedFrame);
	vEvent.push_back(param->hStatsUpdate);
	vEvent.push_back(param->hIn_AOCGeneralTrigger);
	vEvent.push_back(param->hIn_PupMask_available);
	vEvent.push_back(param->hIn_ZER_cal_ready);
	vEvent.push_back(param->hIn_ZON_cal_ready);
	vEvent.push_back(param->hIn_DM_static_map_updated);
	vEvent.push_back(param->hIn_RefCentroidsSet);

	//map<string, HANDLE>::iterator it = param->mapHandle.begin();
	//for (;it!=param->mapHandle.end();++it)
	//{
	//	vEvent.push_back((*it).second);
	//	cout << "Registering '" << (*it).first << "' handle.\n";
	//}
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	while(param->bEndThread == false)
	{
		WaitResult = WaitForMultipleObjects(nEvents,&vEvent[0],FALSE,param->WaitTimeout);

		switch(WaitResult)
		{
			case WAIT_OBJECT_0:
			{
				// TERMINATE THREAD
				param->bEndThread = true;
				break;
			}
			case WAIT_OBJECT_0+1:
			{
				// MISSED FRAME SIGNAL
				cout << "Missed frame.\n";
				break;
			}
			case WAIT_OBJECT_0+2:
			{
				// SIGNAL THE DISPLAY IMAGE BUFFER THAT A NEW IMAGE IS WRITTEN
				SHM_imarray[0].md[0].write = 0;				// signal done writing data
				SHM_imarray[0].md[0].cnt0 ++;				// increment counter!
				SHM_imarray[0].md[0].cnt1 ++;				// STATISTICS UPDATE SIGNAL

				SHM_Slopes[0].md[0].write = 0;				// signal done writing data
				SHM_Slopes[0].md[0].cnt0 ++;					// increment counter!
				SHM_Slopes[0].md[0].cnt1 ++;					// STATISTICS UPDATE SIGNAL

				SHM_DM_state[0].md[0].write = 0;				// signal done writing data
				SHM_DM_state[0].md[0].cnt0 ++;					// increment counter!
				SHM_DM_state[0].md[0].cnt1 ++;					// STATISTICS UPDATE SIGNAL

				cout << string_format("count, time_interval, FPS >> %.3f %d %.3es %.3f %.6f\n",TP_iXon.display_image_interval_s, TP_iXon.display_image_count,TP_iXon.time_diff_s,TP_iXon.FPS, 1 / TP_iXon.FPS);

				break;
			}
			case WAIT_OBJECT_0+3:
			{
				// _AOC_GENERAL_TRIGGER
				cout << "Detected '_AOC_GENERAL_TRIGGER' event!\n";

				break;
			}
			case WAIT_OBJECT_0+4:
			{
				// _AOC_PUPMASK_AVAILABLE
				cout << "Detected '_AOC_PUPMASK_AVAILABLE' event!\n";
				// GET THE SHARED MEMORY FILE DATA
				int *pstart = SHM_DM_slope_mask->array.SI32;
				int *pend = SHM_DM_slope_mask->array.SI32 + SHM_DM_slope_mask->md->nelement;
				TP_calcSH.vSH_PupilMask.resize(SHM_DM_slope_mask->md->nelement);
				//std::copy(pstart, pend, &TP_calcSH.vSH_PupilMask[0]);
				// SET THE MAPPING OF THE VALID SUBAPERTURES

				break;
			}
			case WAIT_OBJECT_0+5:
			{
				// _AOC_ZER_CAL_READY
				cout << "Detected '_AOC_ZER_CAL_READY' event!\n";

				// GET THE DATA THROUGH FITS FILE
				vector<int> vSize;
				string strVal;
				string strComment;

				FITS.OpenFITS("R:\\ZER_CAL.fits");
				FITS.GetKeywordValue("NMODES",strVal,strComment);
				TP_calcSH.nModes_ZER = std::atoi(strVal.data());
				TP_calcSH.vDM_Modes_ZER_Slope.resize(TP_calcSH.nModes_ZER * TP_calcSH.nx*TP_calcSH.ny * 2);
				TP_calcSH.vDM_Modes_ZER_iNorm.resize(TP_calcSH.nModes_ZER);
				TP_calcSH.vDM_ProjOnModes_ZER.resize(TP_calcSH.nModes_ZER);

				FITS.GetKeywordValue("NCALFRM",strVal,strComment);
				TP_calcSH.nCalFrameCount = std::atoi(strVal.data());

				cout << string_format("Calibrating with %d ZERnike modes (each with %d recorded frames).\n", TP_calcSH.nModes_ZER, TP_calcSH.nCalFrameCount);

				int retval = FITS.CloseFITS();

				if (TP_calcSH.nCalFrameCount <= 0)
				{
					// ERROR
					cout << string_format("Invalid calibration frame count per mode (%d).\n",TP_calcSH.nCalFrameCount);
					TP_calcSH.bCalMode = false;
					break;
				}

				// LOAD THE nModes_ZER*nAct ARRAY (DIRECT DM COMMAND FOR EACH MODE, NOT THE nact_lin_x*nact_lin_y SQUARE DM)
				if (FITS.Load("R:\\ZER_CAL.fits",TP_calcSH.vDM_Cmd_ZER,vSize) != SUCCESS)
				{
					// ERROR
					cout << "Failed to load 'R:\\ZER_CAL.fits'.\n";
					TP_calcSH.bCalMode = false;
					break;
				}

				if ((vSize.size() != 2) || (vSize[1] != TP_calcSH.nModes_ZER) || (vSize[0] != TP_calcSH.nAct))
				{
					// ERROR
					cout << string_format("Inconsistent data sizes (FITS: %dx%d, System: %dx%d).\n",
						vSize[0],vSize[1],
						TP_calcSH.nModes_ZER,
						TP_calcSH.nAct);
					TP_calcSH.bCalMode = false;
					break;
				}

				cout << string_format("ZERnike calibration data loaded (nModes: %d).\n",TP_calcSH.nModes_ZER);

				TP_calcSH.bCalDataAvailable_ZER = true;
				TP_calcSH.bCalMode = true;

				// SET THE '_AOC_OK' EVENT TO SIGNAL THAT THE PROCESSING HAS BEEN DONE
				SetEvent(param->hOut_OK);

				break;
			}
			case WAIT_OBJECT_0+6:
			{
				// _AOC_ZON_CAL_READY
				cout << "Detected '_AOC_ZON_CAL_READY' event!\n";

				// GET THE DATA THROUGH FITS FILE
				vector<int> vSize;
				string strVal;
				string strComment;

				FITS.OpenFITS("R:\\ZON_CAL.fits");
				FITS.GetKeywordValue("NMODES",strVal,strComment);
				TP_calcSH.nModes_ZON = std::atoi(strVal.data());
				TP_calcSH.vDM_Modes_ZON_Slope.resize(TP_calcSH.nModes_ZON * TP_calcSH.nx*TP_calcSH.ny * 2);
				TP_calcSH.vDM_Modes_ZON_iNorm.resize(TP_calcSH.nModes_ZON);
				TP_calcSH.vDM_ProjOnModes_ZON.resize(TP_calcSH.nModes_ZON);

				FITS.GetKeywordValue("NCALFRM",strVal,strComment);
				TP_calcSH.nCalFrameCount = std::atoi(strVal.data());

				cout << string_format("Calibrating with %d ZONal modes (each with %d recorded frames).", TP_calcSH.nModes_ZON, TP_calcSH.nCalFrameCount);

				int retval = FITS.CloseFITS();

				if (TP_calcSH.nCalFrameCount <= 0)
				{
					// ERROR
					cout << string_format("Invalid calibration frame count per mode (%d).\n",TP_calcSH.nCalFrameCount);
					TP_calcSH.bCalMode = false;
					break;
				}

				vector<float> vArray;

				// LOAD THE nModes_ZON*nAct ARRAY (DIRECT DM COMMAND FOR EACH MODE, NOT THE nact_lin_x*nact_lin_y SQUARE DM)
				if (FITS.Load("R:\\ZON_CAL.fits",TP_calcSH.vDM_Cmd_ZON,vSize) != SUCCESS)
				{
					// ERROR
					cout << "Failed to load 'R:\\ZON_CAL.fits'.\n";
					TP_calcSH.bCalMode = false;
					break;
				}

				if ((vSize.size() != 2) || (vSize[1] != TP_calcSH.nModes_ZON) || (vSize[0] != TP_calcSH.nAct))
				{
					// ERROR
					cout << string_format("Inconsistent data sizes (FITS: %dx%d, System: %dx%d).\n",
						vSize[0],vSize[1],
						TP_calcSH.nModes_ZON,
						TP_calcSH.nAct);
					TP_calcSH.bCalMode = false;
					break;
				}

				cout << string_format("ZERnike calibration data loaded (nModes: %d).\n",TP_calcSH.nModes_ZON);


				TP_calcSH.bCalDataAvailable_ZON = true;
				TP_calcSH.bCalMode = true;

				break;
			}
			case WAIT_OBJECT_0+7:
			{
				// _AOC_DM_STATIC_MAP_UPDATED
				cout << "Detected '_AOC_DM_STATIC_MAP_UPDATED' event!\n";
				// update_DM_static_buffer DOES *NOT* UPDATE THE DM COMMAND MAP (IT JUST UPDATES THE STATIC MAP)
				//update_DM_static_buffer();
				//collapse_DM_buffers();
				break;
			}
			case WAIT_OBJECT_0+8:
			{
				// _AOC_REF_CENTROIDS_SET
				cout << "Detected '_AOC_REF_CENTROIDS_SET' event!\n";
				ApplyRefCentroids();
				break;
			}
			case WAIT_TIMEOUT:
			{
				// TIMEOUT
				break;
			}
		}
	}

	return 0;
}
