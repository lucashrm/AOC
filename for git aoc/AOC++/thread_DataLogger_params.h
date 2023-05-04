#pragma once

#include <cstdlib>
#include <Windows.h>

#include <MCODE.h>
#include <MCODE_Export.h>
#include <MCODE_Export_FITS.h>

using namespace MultiCODE;

typedef struct
{
	bool			bEndThread;
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;
	HANDLE			hInitializeLogBufSize;
	HANDLE			hStartLog;
	HANDLE			hStopLog;

	// TRIGERING EVENTS
	HANDLE			hDataBlockAvailable;

	size_t			nCircBuf;

	size_t			nBlocksToLog;
	size_t			dataBlockSize;
	size_t			dataLineSize;
	vector<vector<float>>	vLogBufferDataSH;

	size_t			LoggerCounter;
	size_t			circBufferLoggerCounter;
	size_t			circBufferCalcSHCounter;
} ThreadParams_DataLogger;





typedef struct _ThreadParams_DataSaver
{
	bool			bEndThread;
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;
	HANDLE			hSaveData;

	ThreadParams_DataLogger
					*TP_DataLogger;

	string			strLogDir;
	string			strLogBasename;
	string			strLogDate;
	string			strFilename;

	MCODE_Export_FITS FITS;

	_ThreadParams_DataSaver()
	{
		// CONSTRUCTOR
		strLogDir = "D:\\DATA_AOC\\";
		strLogBasename = "AOC_Log";
		strLogDate = "YYYY.MM.DDTHRhMnmSS.mmms";
	};

} ThreadParams_DataSaver;

mcINT16 StartThread_DataLogger();
mcINT16 TerminateThread_DataLogger();
static DWORD Thread_DataLogger(void *pParam);

mcINT16 SetTimeString(string& out);
mcINT16 StartThread_DataSaver();
mcINT16 TerminateThread_DataSaver();
static DWORD Thread_DataSaver(void *pParam);

