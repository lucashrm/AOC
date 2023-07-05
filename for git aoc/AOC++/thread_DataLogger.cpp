#include <cstdlib>
#include <vector>
#include <iostream>

#include <MCODE.h>
#include <MCODE_Export.h>
#include <MCODE_Export_FITS.h>

#include "ImageStruct.h"
#include "thread_CalcSH_params.h"
#include "thread_DataLogger_params.h"
#include <direct.h>


using namespace std;
using namespace MultiCODE;

extern MultiCODE::MCODE_LogManager* Log;
extern string						LogMsg;


extern int						N_DM_COMMAND_MAPS;

extern ThreadParams_calcSH		TP_calcSH;
extern DWORD 					Thread_calcSH_ID;
extern HANDLE					hThread_calcSH;

extern ThreadParams_DataLogger	TP_DataLogger;
extern DWORD 					Thread_DataLogger_ID;
extern HANDLE					hThread_DataLogger;

extern ThreadParams_DataSaver	TP_DataSaver;
extern DWORD 					Thread_DataSaver_ID;
extern HANDLE					hThread_DataSaver;

static int						nThreads = 0;
static HANDLE					hThreadManager = INVALID_HANDLE_VALUE;
static int						ThreadManagerID = 0;

extern vector<float>			vtmpLO;
extern vector<float>			vtmpHO;
extern vector<double>			vTT_drift;
extern vector<double>			vHO_drift;

extern float					leakage_gain;

extern vector<float>			vCmdMatrix;

extern IMAGE*					SHM_imarray;			// SHARED MEMORY INSTANT IXON IMAGE
extern IMAGE*					SHM_Slopes;				// SHARED MEMORY N FRAMES HISTORY DATA (DISPLAY RATE)
extern IMAGE*					SHM_RefCent;			// SHARED MEMORY REFERENCE CENTROIDS DATA
extern IMAGE*					SHM_DM_slope_mask;		// SHARED MEMORY PUPIL MASK DATA
extern IMAGE*					SHM_DM_state;			// SHARED MEMORY DM STATE
extern IMAGE*					SHM_ModesInfo;			// SHARED MEMORY MODE INFORMATION
extern vector<IMAGE*>			vSHM_DM_BUFFERS;		// SHARED MEMORY

extern vector<vector<float>>	vDM_cmd_buffers;		// DM CORRECTION BUFFERS (N_DM_COMMAND_MAPS x nAct)
extern vector<float>			vDM_static_cmd;			// DM STATIC COMMAND MAP (MAP CORRESPONDING TO NON-TIME-CRITICAL MAPS)
extern vector<size_t>			vCollapseStaticInd;		// DM STATIC MAPS INDICES
extern vector<size_t>			vCollapseDynamicInd;	// DM DYNAMIC (TIME-CRITICAL) INDICES = NOT(vCollapseStaticInd)
extern vector<float>			vDM_cur_cmd;			// DM INSTANTANEOUS COMMAND MAP = vDM_cmd_buffers[vCollapseIndices] + vDM_static_cmd
extern vector<long long>		vTimeStamp;				// TIME VECTOR RECORDED JUST BEFORE SETTING THE DM INTO POSITION

extern vector<int>				vDMSquare2DMmap;

// VALID SUBAPERTURE INDEX
extern vector<int>				vSH_ValidCell;


//LARGE_INTEGER t0;
//LARGE_INTEGER t1;
//
//LARGE_INTEGER liFreq;

static DWORD Thread_DataLogger(void* pParam)
{
	ThreadParams_DataLogger* param = (ThreadParams_DataLogger*)pParam;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hInitializeLogBufSize);
	vEvent.push_back(param->hStopLog);
	vEvent.push_back(param->hStartLog);
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	LARGE_INTEGER t0;
	LARGE_INTEGER now;

	LARGE_INTEGER liCPUFreq;
	QueryPerformanceFrequency(&liCPUFreq);

	double tick_period_microsecond = 1. / double(liCPUFreq.QuadPart) * 1.e6;

	int waitResult;

	while (param->bEndThread == false)
	{
		waitResult = WaitForMultipleObjects(nEvents, &vEvent[0], FALSE, INFINITE);

		switch (waitResult)
		{
		case WAIT_OBJECT_0:
		{
			param->bEndThread = true;
			break;
		}
		case WAIT_OBJECT_0 + 1:
		{
			// INITIALIZE LOG BUFFERS
			cout << string_format("DataLogger log buffer size initialization.\n");
			vector<float> vTmp(param->dataLineSize*param->dataBlockSize*param->nBlocksToLog);
			param->vLogBufferDataSH.resize(param->nCircBuf, vTmp);
			break;
		}
		case WAIT_OBJECT_0 + 2:
		{
			// STOP LOGGING PROCESS
			cout << string_format("DataLogger thread successfully stopped.\n");
			ResetEvent(param->hStartLog);
			ResetEvent(param->hStopLog);
			TP_DataLogger.LoggerCounter = 0;
			break;
		}
		case WAIT_OBJECT_0 + 3:
		{
			size_t logSize = param->dataLineSize*param->dataBlockSize*param->nBlocksToLog;
			size_t blockSize = param->dataLineSize*param->dataBlockSize;
			vector<float> vTmp(param->dataLineSize*param->dataBlockSize*param->nBlocksToLog);
			vector<double> vTmpL(param->dataLineSize * param->dataBlockSize * param->nBlocksToLog);
			param->vLogBufferDataSH.resize(param->nCircBuf, vTmp);
			param->vLogBufferDMcmd.resize(param->nCircBuf, vTmp);
			param->vLogBufferTimestamps.resize(param->nCircBuf, vTmpL);

			TP_DataLogger.LoggerCounter = 0;
			size_t blockCounter = 0;
			bool bStop = false;

			vector<HANDLE> vEvent;
			vEvent.push_back(param->hStopLog);
			vEvent.push_back(param->hDataBlockAvailable);

			// START LOGGING PROCESS

			while (true)
			{
				int waitResult = WaitForMultipleObjects(vEvent.size(),&vEvent[0],FALSE,INFINITE);
				
				switch(waitResult)
				{
				case WAIT_OBJECT_0:
					bStop = true;
					break;
				case WAIT_OBJECT_0+1:
					std::copy(TP_calcSH.vCircBuf_FlatSHData[param->circBufferCalcSHCounter].begin(), TP_calcSH.vCircBuf_FlatSHData[param->circBufferCalcSHCounter].end(),
						&param->vLogBufferDataSH[param->circBufferLoggerCounter][blockCounter*blockSize]);
					std::copy(TP_calcSH.vCircBuf_DM_cmd[param->circBufferCalcSHCounter].begin(), TP_calcSH.vCircBuf_DM_cmd[param->circBufferCalcSHCounter].end(),
						&param->vLogBufferDMcmd[param->circBufferLoggerCounter][blockCounter * blockSize]);
					std::copy(TP_calcSH.vCircBuf_Timestamps[param->circBufferCalcSHCounter].begin(), TP_calcSH.vCircBuf_Timestamps[param->circBufferCalcSHCounter].end(),
						&param->vLogBufferTimestamps[param->circBufferLoggerCounter][blockCounter * blockSize]);
					break;
				}

				if (bStop == true)
				{
					break;
				}

				if (++blockCounter >= param->nBlocksToLog)
				{
					++param->LoggerCounter;
					param->circBufferLoggerCounter = param->LoggerCounter % param->nCircBuf;
					blockCounter = 0;
					SetEvent(TP_DataSaver.hSaveData);
				}
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

	CloseHandle(TP_DataLogger.hTerminateThread); TP_DataLogger.hTerminateThread = INVALID_HANDLE_VALUE;
	CloseHandle(TP_DataLogger.hStartLog); TP_DataLogger.hStartLog = INVALID_HANDLE_VALUE;
	CloseHandle(TP_DataLogger.hStopLog); TP_DataLogger.hStopLog = INVALID_HANDLE_VALUE;
	CloseHandle(hThread_DataLogger); hThread_DataLogger = INVALID_HANDLE_VALUE;

	return 0;
}

mcINT16 StartThread_DataLogger()
{
	string threadName = "Thread_DataLogger";

	TP_DataLogger.hThreadCreated = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_DataLogger.hTerminateThread = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_DataLogger.hInitializeLogBufSize = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_DataLogger.hStartLog = CreateEvent(NULL, TRUE, FALSE, NULL);
	TP_DataLogger.hStopLog = CreateEvent(NULL, TRUE, FALSE, NULL);
	TP_DataLogger.hDataBlockAvailable = CreateEvent(NULL, FALSE, FALSE, NULL);

	TP_DataLogger.nCircBuf = 2;
	TP_DataLogger.dataLineSize = 500;
	TP_DataLogger.dataBlockSize = 50;
	TP_DataLogger.nBlocksToLog = 1000;

	hThread_DataLogger = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)Thread_DataLogger, &TP_DataLogger, THREAD_PRIORITY_NORMAL, &Thread_DataLogger_ID);

	DWORD WaitResult = WaitForSingleObject(TP_DataLogger.hThreadCreated, THREAD_CREATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		LogMsg = string_format("[%s] thread %s created successuflly.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SUCCESS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		break;
	case WAIT_TIMEOUT:
		LogMsg = string_format("[%s] %s thread creation TIMEOUT.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_ABANDONED:
		LogMsg = string_format("[%s] %s thread creation ABANDONED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_FAILED:
		LogMsg = string_format("[%s] %s thread creation FAILED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	}

	return SUCCESS;
}

mcINT16 TerminateThread_DataLogger()
{
	string threadName = "EvtWatch";

	SetEvent(TP_DataLogger.hTerminateThread);

	DWORD WaitResult = WaitForSingleObject(hThread_DataLogger, THREAD_TERMINATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		LogMsg = string_format("[%s] thread %s terminated successuflly.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SUCCESS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		break;
	case WAIT_TIMEOUT:
		LogMsg = string_format("[%s] %s thread termination TIMEOUT.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_ABANDONED:
		LogMsg = string_format("[%s] %s thread termination ABANDONED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_FAILED:
		LogMsg = string_format("[%s] %s thread termination FAILED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	}

	return SUCCESS;
}










mcINT16 SetTimeString(SYSTEMTIME& myST, string& out)
{
	out = string_format("%04d.%02d.%02dT%02dh%02dm%02d.%03ds", myST.wYear, myST.wMonth, myST.wDay,
		myST.wHour, myST.wMinute, myST.wSecond, myST.wMilliseconds);

	return SUCCESS;
}


static DWORD Thread_DataSaver(void* pParam)
{
	ThreadParams_DataSaver* param = (ThreadParams_DataSaver*)pParam;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hSaveData);
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	LARGE_INTEGER t0;
	LARGE_INTEGER now;

	LARGE_INTEGER liCPUFreq;
	QueryPerformanceFrequency(&liCPUFreq);

	double tick_period_microsecond = 1. / double(liCPUFreq.QuadPart) * 1.e6;

	int waitResult;

	while (param->bEndThread == false)
	{
		waitResult = WaitForMultipleObjects(nEvents, &vEvent[0], FALSE, INFINITE);

		switch (waitResult)
		{
		case WAIT_OBJECT_0:
		{
			param->bEndThread = true;
			break;
		}
		case WAIT_OBJECT_0 + 1:
		{
			// SAVE LOG FILE
			SYSTEMTIME myST;
			GetSystemTime(&myST);
			SetTimeString(myST, param->strLogDate);
			
			/*if (!TP_calcSH.bCalMode) {
				std::string strDirName = string_format("%s_%s_%s", "prefix_", param->strLogDate.data(), "_NONE");
				param->strLogDir = "C:\\Users\\lucas\\Documents\\STAGE\\Misc\\log\\";
				_mkdir(param->strLogDir.c_str());
			}*/

			param->strFilename = string_format("%s_%s_%04d", param->strLogDate.data(), param->strLogBasename.data(),
				param->TP_DataLogger->LoggerCounter);

			std::string strFilenameT = string_format("%s_%s_%04d", param->strLogDate.data(), "_timestamp",
				param->TP_DataLogger->LoggerCounter);

			std::string strFilenameDM = string_format("%s_%s_%04d", param->strLogDate.data(), "_dm_cmd",
				param->TP_DataLogger->LoggerCounter);

			cout << string_format("DataSaver saved file '%s' to folder '%s'.\n", param->strFilename.data(), param->strLogDir.data());
			

			ThreadParams_DataLogger& TPDL = *param->TP_DataLogger;

			if (param->FITS.Save(&TPDL.vLogBufferDataSH[TPDL.circBufferLoggerCounter][0], TPDL.dataLineSize, TPDL.dataBlockSize*TPDL.nBlocksToLog, param->strLogDir, param->strFilename) != 0)
				cout << "Couldn't log data sh" << endl;
			if (param->FITS.Save(&TPDL.vLogBufferDMcmd[TPDL.circBufferLoggerCounter][0], TPDL.dataLineSize, TPDL.dataBlockSize * TPDL.nBlocksToLog, param->strLogDir, strFilenameDM) != 0)
				cout << "Couldn't log dm cmd" << endl;
			if (param->FITS.Save(&TPDL.vLogBufferTimestamps[TPDL.circBufferLoggerCounter][0], TPDL.dataLineSize, TPDL.dataBlockSize * TPDL.nBlocksToLog, param->strLogDir, strFilenameT) != 0)
				cout << "Couldn't log timestamps" << endl;
			
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

	CloseHandle(TP_DataSaver.hTerminateThread); TP_DataSaver.hTerminateThread = INVALID_HANDLE_VALUE;
	CloseHandle(TP_DataSaver.hSaveData); TP_DataSaver.hSaveData = INVALID_HANDLE_VALUE;
	CloseHandle(hThread_DataSaver); hThread_DataSaver = INVALID_HANDLE_VALUE;

	return 0;
}

mcINT16 StartThread_DataSaver()
{
	string threadName = "Thread_DataSaver";

	TP_DataSaver.hThreadCreated = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_DataSaver.hTerminateThread = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_DataSaver.hSaveData = CreateEvent(NULL, FALSE, FALSE, NULL);

	TP_DataSaver.strPrefix = "obs";

	TP_DataSaver.TP_DataLogger = &TP_DataLogger;

	/*std::string timeDir;
	SYSTEMTIME myST;
	GetSystemTime(&myST);
	SetTimeString(myST, timeDir);

	if (!TP_calcSH.bCalMode) {
		std::string strDirName = string_format("%s_%s_%s", "prefix_", timeDir.data(), "_NONE");
		TP_DataSaver.strLogDir = "C:\\Users\\lucas\\Documents\\STAGE\\Misc\\log\\" + strDirName + "\\";
		_mkdir(TP_DataSaver.strLogDir.c_str());
	}*/

	hThread_DataSaver = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)Thread_DataSaver, &TP_DataSaver, THREAD_PRIORITY_NORMAL, &Thread_DataSaver_ID);

	DWORD WaitResult = WaitForSingleObject(TP_DataSaver.hThreadCreated, THREAD_CREATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		LogMsg = string_format("[%s] thread %s created successuflly.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SUCCESS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		break;
	case WAIT_TIMEOUT:
		LogMsg = string_format("[%s] %s thread creation TIMEOUT.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_ABANDONED:
		LogMsg = string_format("[%s] %s thread creation ABANDONED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_FAILED:
		LogMsg = string_format("[%s] %s thread creation FAILED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	}

	return SUCCESS;
}

mcINT16 TerminateThread_DataSaver()
{
	string threadName = "DataSaver";

	SetEvent(TP_DataSaver.hTerminateThread);

	DWORD WaitResult = WaitForSingleObject(hThread_DataSaver, THREAD_TERMINATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		LogMsg = string_format("[%s] thread %s terminated successuflly.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SUCCESS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		break;
	case WAIT_TIMEOUT:
		LogMsg = string_format("[%s] %s thread termination TIMEOUT.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_ABANDONED:
		LogMsg = string_format("[%s] %s thread termination ABANDONED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	case WAIT_FAILED:
		LogMsg = string_format("[%s] %s thread termination FAILED.", __FUNCTION__, threadName.data());
		Log->StackAddMsg(LOG_SERIOUS, 0, __LINE__, __FILEREF__, LogMsg);
		cout << LogMsg << endl;
		return SERIOUS;
		break;
	}

	return SUCCESS;
}

