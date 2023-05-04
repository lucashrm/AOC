#include <cstdlib>
#include <vector>
#include <iostream>

#include <MCODE.h>
#include <MCODE_SystemUtilities.h>
#include <MCODE_StdFunc.h>

#include "ImageStruct.h"

#include "thread_DMcmd_params.h"

#include "thread_calcSH_params.h"

#undef SUCCESS
#include <asdkDM.h>

#define SUCCESS unsigned short(0)

using namespace std;
using namespace MultiCODE;

extern MultiCODE::MCODE_LogManager	*Log;
extern string						LogMsg;

extern int					N_DM_COMMAND_MAPS;

extern ThreadParams_DMcmd	TP_DMcmd;
extern DWORD 				Thread_DMcmd_ID;
extern HANDLE				hThread_DMcmd;

static int					nThreads = 0;
static HANDLE				hThreadManager = INVALID_HANDLE_VALUE;
static int					ThreadManagerID = 0;

extern IMAGE				*SHM_imarray;			// SHARED MEMORY INSTANT IXON IMAGE
extern IMAGE				*SHM_Slopes;			// SHARED MEMORY N FRAMES HISTORY DATA (DISPLAY RATE)
extern IMAGE				*SHM_RefCent;			// SHARED MEMORY REFERENCE CENTROIDS DATA
extern IMAGE				*SHM_DM_slope_mask;		// SHARED MEMORY PUPIL MASK DATA
extern IMAGE				*SHM_DM_state;			// SHARED MEMORY DM STATE
extern IMAGE				*SHM_ModesInfo;			// SHARED MEMORY MODE INFORMATION
extern vector<IMAGE*>		vSHM_DM_BUFFERS;		// SHARED MEMORY

extern vector<vector<float>>	vDM_cmd_buffers;	// DM CORRECTION BUFFERS (N_DM_COMMAND_MAPS x nAct)
extern vector<float>			vDM_static_cmd;		// DM STATIC COMMAND MAP (MAP CORRESPONDING TO NON-TIME-CRITICAL MAPS)
extern vector<size_t>			vCollapseStaticInd;	// DM STATIC MAPS INDICES
extern vector<size_t>			vCollapseDynamicInd;	// DM DYNAMIC (TIME-CRITICAL) INDICES = NOT(vCollapseStaticInd)
extern vector<float>			vDM_cur_cmd;		// DM INSTANTANEOUS COMMAND MAP = vDM_cmd_buffers[vCollapseIndices] + vDM_static_cmd

extern vector<int>				vDMSquare2DMmap;

extern acs::DM					dm;
extern vector<acs::Scalar>		vDM_cmd;

extern ThreadParams_calcSH	TP_calcSH;
extern DWORD 				Thread_calcSH_ID;
extern HANDLE				hThread_calcSH;

int collapse_DM_buffers()
{
	vector<float>::iterator itA = vDM_cur_cmd.begin();

	for (size_t k=0;k<TP_calcSH.nAct;++k, ++itA)
	{
		*itA = vDM_static_cmd[k];
		for (size_t u=0;u<vCollapseDynamicInd.size();++u)
		{
			float value = *(vSHM_DM_BUFFERS[vCollapseDynamicInd[u]]->array.F + k);
			vDM_cmd_buffers[vCollapseDynamicInd[u]][k] = value;
			*itA += value;
		}
	}

	// SET MAP TO THE LAST SHARED MEMORY BUFFER
	std::copy(vDM_cur_cmd.begin(),vDM_cur_cmd.end(),SHM_DM_state->array.F);
	SHM_DM_state->md[0].cnt0++;

	return 0;
}

int update_DM_static_buffer()
{
	for (size_t k=0;k<vCollapseStaticInd.size();++k)
	{
		float *pstart = vSHM_DM_BUFFERS[vCollapseStaticInd[k]]->array.F;
		float *pend = vSHM_DM_BUFFERS[vCollapseStaticInd[k]]->array.F + vSHM_DM_BUFFERS[vCollapseStaticInd[k]]->md->nelement;
		// WARNING: SHM OBJECTS HAVE A nAct SIZE (NOT nact_lin_x*nact_lin_y)

		std::copy(pstart, pend, vDM_cmd_buffers[vCollapseStaticInd[k]].begin());
	}

	vector<float>::iterator itA = vDM_static_cmd.begin();

	for (size_t k=0;k<TP_calcSH.nAct;++k,++itA)
	{
		*itA = vDM_cmd_buffers[vCollapseStaticInd[0]][k];
		for (size_t u=1;u<vCollapseStaticInd.size();++u)
		{
			*itA += vDM_cmd_buffers[vCollapseStaticInd[u]][k];
		}
	}

	return 0;
}

int DM_set_to_position()
{
	vector<float>::iterator itf = vDM_cur_cmd.begin();
	vector<double>::iterator itd = vDM_cmd.begin();
	for (size_t k = 0; k < vDM_cmd.size(); ++k, ++itf, ++itd)
	{
		*itd = double(*itf);
	}

	if (TP_DMcmd.bEnableDM == true)
	{
		try
		{
			//dm.Reset();
			dm.Send(vDM_cmd.data());
		}
		catch (std::exception e)
		{
			cout << "'DM_set_to_position' EXCEPTION: " << e.what() << endl;
			return -1;
		}
	}

	return 0;
}

static DWORD Thread_DMcmd(void *pParam)
{
	ThreadParams_DMcmd *param = (ThreadParams_DMcmd*)pParam;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hStopServer);
	vEvent.push_back(param->hStartServer);
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	LARGE_INTEGER t0;
	LARGE_INTEGER now;

	LARGE_INTEGER liCPUFreq;
	QueryPerformanceFrequency(&liCPUFreq);

	double tick_period_microsecond = 1./double(liCPUFreq.QuadPart) * 1.e6;

	int waitResult;

	while(param->bEndThread == false)
	{
		waitResult = WaitForMultipleObjects(nEvents,&vEvent[0],FALSE,param->waitTimeout_ms);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			{
				param->bEndThread = true;
				break;
			}
		case WAIT_OBJECT_0+1:
			{
				// STOP DM SERVER
				param->waitTimeout_ms = INFINITE;
				ResetEvent(param->hStopServer);
				cout << string_format("DMcmd thread successfully suspended.");
				break;
			}
		case WAIT_OBJECT_0+2:
			{
				// START DM SERVER
				QueryPerformanceCounter(&t0);
				collapse_DM_buffers();
				DM_set_to_position();

				for (;;)
				{
					if (WaitForSingleObject(param->hStopServer, 0) == WAIT_OBJECT_0)
					{
						cout << "'hStopServer' event was triggered. Stopping loop.\n";
						break;
					}

					

					QueryPerformanceCounter(&now);

					int count = 0;
					while (((double(now.QuadPart) - double(t0.QuadPart)) * tick_period_microsecond) < param->loopTiming_microsec)
					{
						collapse_DM_buffers();
						++count;
						QueryPerformanceCounter(&now);
					}

					DM_set_to_position();

					QueryPerformanceCounter(&t0);
				}

				break;
			}
		case WAIT_FAILED:
			{
				cout << string_format("WAIT_FAILED!!\n");
				break;
			}
		case WAIT_ABANDONED:
			{
				cout << string_format("WAIT_ABANDONED!!\n");
				break;
			}
		}
	}

	CloseHandle(TP_DMcmd.hTerminateThread); TP_DMcmd.hTerminateThread = INVALID_HANDLE_VALUE;
	CloseHandle(TP_DMcmd.hStartServer); TP_DMcmd.hStartServer = INVALID_HANDLE_VALUE;
	CloseHandle(TP_DMcmd.hStopServer); TP_DMcmd.hStopServer = INVALID_HANDLE_VALUE;
	CloseHandle(hThread_DMcmd); hThread_DMcmd = INVALID_HANDLE_VALUE;

	return 0;
}

mcINT16 StartThread_DMcmd()
{
	string threadName = "Thread_DMcmd";

	TP_DMcmd.hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_DMcmd.hTerminateThread = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_DMcmd.hStartServer = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_DMcmd.hStopServer = CreateEvent(NULL,TRUE,FALSE,NULL);

	TP_DMcmd.bEndThread = false;

	TP_DMcmd.bEnableDM = false;

	TP_DMcmd.waitTimeout_ms = INFINITE;

	TP_DMcmd.loopTiming_microsec = 700;

	hThread_DMcmd = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_DMcmd,&TP_DMcmd,THREAD_PRIORITY_HIGHEST,&Thread_DMcmd_ID);

	cout << string_format("HANDLE for thread: %x  %ud\n",hThread_DMcmd,Thread_DMcmd_ID);

	DWORD WaitResult = WaitForSingleObject(TP_DMcmd.hThreadCreated,THREAD_CREATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		CloseHandle(TP_DMcmd.hThreadCreated);
		TP_DMcmd.hThreadCreated = INVALID_HANDLE_VALUE;
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

mcINT16 TerminateThread_DMcmd()
{
	string threadName = "DMcmd";

	SetEvent(TP_DMcmd.hStopServer);
	SetEvent(TP_DMcmd.hTerminateThread);

	DWORD WaitResult = WaitForSingleObject(hThread_DMcmd,THREAD_TERMINATION_TIMEOUT);

	switch(WaitResult)
	{
	case WAIT_OBJECT_0:
		hThread_DMcmd = INVALID_HANDLE_VALUE;
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
