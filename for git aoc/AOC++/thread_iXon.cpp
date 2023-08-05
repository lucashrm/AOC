#include <iostream>
#include <cstdlib>

#include <MultiCODE.h>
#include <MCODE.h>

#include <MCODE_SystemUtilities.h>
#include <MCODE_StdFunc.h>

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

#include "thread_iXon_params.h"
#include "thread_calcSH_params.h"
#include "thread_EvtWatch_params.h"

extern MultiCODE::MCODE_LogManager	*Log;
extern string						LogMsg;

extern ThreadParams_iXon			TP_iXon;
extern ThreadParams_Simu_iXon		TP_Simu_iXon;
extern DWORD						Thread_iXon_ID;
extern DWORD						Thread_Simu_iXon_ID;
extern HANDLE						hThread_iXon;
extern HANDLE						hThread_Simu_iXon;

extern ThreadParams_calcSH			TP_calcSH;
extern DWORD						Thread_calcSH_ID;
extern HANDLE						hThread_calcSH;

extern ThreadParams_EvtWatch		TP_EvtWatch;
extern DWORD						Thread_EvtWatch_ID;
extern HANDLE						hThread_EvtWatch;

extern IMAGE						*SHM_imarray;		// pointer to image

extern int							localCount;

#pragma comment(lib,"atmcd64m.lib")

using namespace std;

mcINT16 StartThread_iXon()
{
	string threadName = "iXon";

	// THREAD CREATION/TERMINATION EVENTS
	TP_iXon.hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_iXon.hTerminateThread = CreateEvent(NULL,FALSE,FALSE,NULL);

	// TRIGGERING EVENTS
	TP_iXon.hStartAcquisition = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_iXon.hStopAcquisition = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_iXon.hShutdownCamera = CreateEvent(NULL,FALSE,FALSE,NULL);
	
	// SIGNALING EVENTS
	TP_iXon.hFrameReady = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_iXon.hAcquisitionDone = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_iXon.hLiveStream = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_iXon.hExt_StatsUpdate = TP_EvtWatch.hStatsUpdate;
	TP_iXon.hExt_MissedFrame = TP_EvtWatch.hMissedFrame;
	TP_iXon.hGetLastImage = CreateEvent(NULL, FALSE, FALSE, NULL);

	// THREAD VARIABLE INITIALIZATION
	//TP_iXon.FPSinterval = 1000;
	TP_iXon.display_image_interval_s = 0.1;
	TP_iXon.display_image_count = 50;//int(1. / TP_iXon.kinTime_s * TP_iXon.display_image_interval_s + 0.5);
	TP_iXon.bEndThread = false;
	TP_iXon.bStreamON = false;
	TP_iXon.acqTimeout_ms = 1000;
	TP_iXon.WaitTimeout = INFINITE;
	TP_iXon.bCamInitialized = false;
	TP_iXon.camID = -1;

	TP_iXon.bAsyncProc = false;

	TP_iXon.nCalFrames = 0;
	TP_iXon.bCalibration = false;

	TP_iXon.testLocalCount = 0;

	TP_iXon.vImgBuffer.resize(TP_calcSH.DX * TP_calcSH.DY);

	hThread_iXon = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_iXon,&TP_iXon,THREAD_PRIORITY_ABOVE_NORMAL,&Thread_iXon_ID);

	DWORD WaitResult = WaitForSingleObject(TP_iXon.hThreadCreated,THREAD_CREATION_TIMEOUT);

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

mcINT16 TerminateThread_iXon()
{
	string threadName = "iXon";

	SetEvent(TP_iXon.hTerminateThread);

	DWORD WaitResult = WaitForSingleObject(hThread_iXon,THREAD_TERMINATION_TIMEOUT);

	switch(WaitResult)
	{
	case WAIT_OBJECT_0:
		hThread_iXon = INVALID_HANDLE_VALUE;
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

// ==========================================================================
//                         SIMULATION IXON THREAD
// ==========================================================================


mcINT16 StartThread_Simu_iXon()
{
	string threadName = "Simu_iXon";

	DWORD waitResult;

	// THREAD CREATION/TERMINATION EVENTS
	TP_Simu_iXon.hThreadCreated = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_Simu_iXon.hTerminateThread = CreateEvent(NULL, FALSE, FALSE, NULL);

	// THREAD TRIGGER EVENTS
	TP_Simu_iXon.hStartLiveStream = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_Simu_iXon.hStopLiveStream = CreateEvent(NULL, FALSE, FALSE, NULL);

	// SIGNALING EVENTS
	TP_Simu_iXon.hLiveStreamDone = CreateEvent(NULL, FALSE, FALSE, NULL);
	TP_Simu_iXon.hMissedFrame = CreateEvent(NULL, FALSE, FALSE, NULL);

	// PARAMETERS
	TP_Simu_iXon.bEndThread = false;
	TP_Simu_iXon.WaitTimeout = INFINITE;
	TP_Simu_iXon.currentIndex = 0;
	TP_Simu_iXon.fitsDir = "C:\\Users\\lucas\\Documents\\STAGE\\FITS\\simu_square_2\\";

	hThread_Simu_iXon = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)Thread_Simu_iXon, &TP_Simu_iXon, THREAD_PRIORITY_ABOVE_NORMAL, &Thread_Simu_iXon_ID);

	waitResult = WaitForSingleObject(TP_Simu_iXon.hThreadCreated, THREAD_CREATION_TIMEOUT);

	switch (waitResult)
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

mcINT16 TerminateThread_Simu_iXon()
{
	string threadName = "Simu_iXon";

	SetEvent(TP_Simu_iXon.hTerminateThread);

	DWORD WaitResult = WaitForSingleObject(hThread_Simu_iXon, THREAD_TERMINATION_TIMEOUT);

	switch (WaitResult)
	{
	case WAIT_OBJECT_0:
		hThread_iXon = INVALID_HANDLE_VALUE;
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

void Simulation(ThreadParams_Simu_iXon* param)
{
	SHM_imarray[0].md[0].write = 1;
	std::copy(SimulationLoader::getInstance().get(param->currentIndex).begin(),
		SimulationLoader::getInstance().get(param->currentIndex).end(), TP_iXon.vImgBuffer.begin());
	if (SimulationLoader::getInstance().get(param->currentIndex).size() == 0)
		cout << "non" << endl;
	//cout << param->currentIndex << endl;
	//MCODE_Export_FITS fits;
	//fits.Save(SHM_imarray[0].array.SI32, 128, 128, "C:\\Users\\lucas\\Documents\\STAGE\\FITS\\", "test");
	param->currentIndex++;
	if (param->currentIndex >= param->maxImages)
		param->currentIndex = 0;
	param->imgCount++;
	//SHM_imarray[0].md[0].cnt0++;
}

mcUINT32 Thread_Simu_iXon(void *pParam)
{
	ThreadParams_Simu_iXon* param = (ThreadParams_Simu_iXon*)pParam;

	DWORD waitResult;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hStartLiveStream);
	vEvent.push_back(param->hStopLiveStream);
	vEvent.push_back(param->hMissedFrame);

	SetEvent(param->hThreadCreated);

	SimulationLoader::getInstance().loadFullRepo(param->fitsDir);

	param->maxImages = SimulationLoader::getInstance().getAllValues().size();
	param->currentIndex = 0;
	param->imgCount = 0;

	LARGE_INTEGER li_freq;
	LARGE_INTEGER li_t0;
	LARGE_INTEGER li_t1;
	LARGE_INTEGER li_prev;
	//LARGE_INTEGER li_now;
	double pc_freq = 0.;
	double pc_start = 0.;

	TP_Simu_iXon.isMissed = false;

	QueryPerformanceFrequency(&li_freq);

	double tick_period_microsecond = 1. / double(li_freq.QuadPart) * 1.e6;

	double dt = 0;

	while (!param->bEndThread) {
		waitResult = WaitForMultipleObjects(vEvent.size(), &vEvent[0], FALSE, param->WaitTimeout);

		switch (waitResult) {
		case WAIT_OBJECT_0:
		{
			param->bEndThread = true;
			break;
		}
		case WAIT_OBJECT_0 + 1:
		{
			param->WaitTimeout = 0;
			param->imgCount = 0;
			TP_iXon.DX = TP_calcSH.DX;
			TP_iXon.DY = TP_calcSH.DY;
			QueryPerformanceCounter(&li_t0);
			break;
		}
		case WAIT_OBJECT_0 + 2:
		{
			param->WaitTimeout = INFINITE;
			SetEvent(param->hLiveStreamDone);
			break;
		}
		case WAIT_OBJECT_0 + 3:
		{
			TP_Simu_iXon.isMissed = true;
			// cout << "is missed in simu" << endl;
		}
		case WAIT_TIMEOUT:
		{
			if (!TP_Simu_iXon.isMissed) {
				QueryPerformanceCounter(&li_t1);
				dt = double((li_t1.QuadPart - li_t0.QuadPart)) * tick_period_microsecond;
				while (dt <= TP_iXon.expTime_s * 1000000) {
					QueryPerformanceCounter(&li_t1);
					dt = double((li_t1.QuadPart - li_t0.QuadPart)) * tick_period_microsecond;
				}
				dt = 0;
				li_t0 = li_t1;
				Simulation(param);
				SetEvent(TP_iXon.hGetLastImage);
			}
			else {
				cout << "is missed in simu" << endl;
			}
			TP_Simu_iXon.isMissed = false;
			break;
		}
		}
	}
	return SUCCESS;
}

mcINT16 CallWaitForAcquisitionTimeout(ThreadParams_iXon *param)
{
	if (!param->bSimulation)
	{
		if (WaitForAcquisitionTimeOut(param->acqTimeout_ms) == DRV_SUCCESS)
			return SUCCESS;
	}
	else
	{
		return SUCCESS;
	}
}

mcINT16 CallGetLastNumberImages(ThreadParams_Simu_iXon* param, ThreadParams_iXon *parami, long &first_im, long &last_im)
{
	if (!parami->bSimulation)
	{
		GetNumberNewImages(&first_im, &last_im);
		return SUCCESS;
	}
	else
	{
		last_im = param->imgCount;
	}
}

mcINT16 GetSimuLastImage(ThreadParams_iXon* param)
{
	DWORD timeout = param->expTime_s < 0.0005 ? 1 : DWORD(2 * param->expTime_s * 1000);
	if (WaitForSingleObject(param->hGetLastImage, timeout) == WAIT_OBJECT_0)
		return SUCCESS;
	else
		return SERIOUS;
}

// ==========================================================================
//                         END SIMULATION IXON THREAD
// ==========================================================================

// ==========================================================================
//                         ACQUISITION THREAD
// ==========================================================================
void* AcquireLoop(ThreadParams_iXon *param)
{
	int status = 0;
	int imgSize = param->DX * param->DY;
	//int nfr = 100, ifr = 0;     // variables used to estimate frame rate
	//float t0 = 0.0, t1 = 0.0;  // time variables for frame rate
	//int idisp = 0, ndisp = 100; // control the text output refresh rate
	double ifrate = 0.;              // frame rate (average over nfr frames)
	double itiming = 0.;
	int local_count = 0;
	localCount = 0;
	int full_count = 0;
	//struct mytimespec now;       // clock readout
	//struct tm* ptm;
	//float *timing = (float*) malloc(nfr * sizeof(float));

	long first_im;				// indices for avail. images in circ buffer
	long last_im;
	long cnt_first;
	long last_read_im;            // index of last read image
	//long first_valid, last_valid; // indices of valid images
	unsigned long error;         // error code

	last_read_im = -1;           // keep track of last read image index
	//for (ifr = 0; ifr < nfr; ifr++) timing[ifr] = 0.1; // init timing array
	//ifr = 1;

	// acquisition loop
	printf("\n");

	if (param->bStreamON)
	{
		::StartAcquisition();
	}
 
	LARGE_INTEGER li_freq;
	LARGE_INTEGER li_t0;
	LARGE_INTEGER li_t1;
	LARGE_INTEGER li_t2;
	LARGE_INTEGER li_prev;
	//LARGE_INTEGER li_now;
	double pc_freq = 0.;
	double pc_start = 0.;

	QueryPerformanceFrequency(&li_freq);
	double iFreq = 1./double(li_freq.QuadPart);
	double avg_delta = 0.;
	double meas_frate = 0.;
	int counter = 0;

	mcINT16 retVal = SUCCESS;

	long delta = 0;
	DWORD camStatus = 0;

	register double time_diff = 0.;

	param->vImgBuffer.resize(imgSize);

	// TAKE A FIRST IMAGE TO RECORD THE COUNTER
	GetNumberNewImages(&first_im, &last_im);
	cnt_first = last_im;

	QueryPerformanceCounter(&param->li_timestamp);
	li_t0 = param->li_timestamp;
	li_prev = param->li_timestamp;

	if (param->bCalibration == true)
	{
		//while (camconf->nleft > 0)
		//while(WaitForAcquisitionTimeOut(param->acqTimeout_ms) == DRV_SUCCESS)
		while(CallWaitForAcquisitionTimeout(param) == SUCCESS)
		{
			QueryPerformanceCounter(&li_t1);
			// CODE FOR 'Internal trigger mode' MODE
		//	error = GetNumberNewImages(&first_im, &last_im);
			CallGetLastNumberImages(&TP_Simu_iXon, param, first_im, last_im);
			//// CODE FOR 'Software trigger' MODE
			//SendSoftwareTrigger();
			//WaitForAcquisition();
			//GetNumberNewImages(&first_im, &last_im);

			delta = last_im - last_read_im;
			if (delta > 1)
			{
				//printf("delta image: %d %d %d %d\n",last_im - last_read_im, last_im, last_read_im,error);
				SetEvent(param->hExt_MissedFrame);
			}
			last_read_im = last_im;
		
			

			if (param->bSimulation) {
				GetSimuLastImage(param);
			}
			else {
				error = GetMostRecentImage(&param->vImgBuffer[0], imgSize);
			}

			if (param->bAsyncProc == true)
			{
				if (TryEnterCriticalSection(&param->CSCopyImgBuffer) == TRUE)
				{
					std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(),param->vImgBufferCopy.begin());
				}
			}
			else
			{
				// FOR TEST PURPOSE!!!!
				std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(),param->vImgBufferCopy.begin());

				// SHARED MEMORY DISPLAY IMAGE BUFFER
				SHM_imarray[0].md[0].write = 1;				// signal you are about to write
				std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(),(long*)SHM_imarray[0].array.SI32);

				retVal = AOC_ProcessSHData(param);
				QueryPerformanceCounter(&li_t2);

				time_diff += double(li_t2.QuadPart - li_t1.QuadPart);

				if (retVal == SUCCESS)
				{

				}
			}

			//SetEvent(param->hFrameReady);

			++full_count;
			++local_count;
			localCount++;


			// CHECK FOR LOOP TERMINATION
			if ((full_count >= param->nCalFrames) || (WaitForSingleObject(param->hStopAcquisition,0) == WAIT_OBJECT_0))
			{
				QueryPerformanceCounter(&li_t2);

				//camconf->nleft = 0;
				if (::AbortAcquisition() != DRV_SUCCESS)
				{
					cout << string_format("'::AbortAcquisition()' FAILED!\n");
				}
				else
				{
					cout << string_format("Acquisition stopped successfully\n\n");
				}
				param->bStreamON = false;

				li_prev = li_t2;
				param->FPS = double(full_count) / (double(li_t2.QuadPart - li_t0.QuadPart) * iFreq);
		
				param->frame_count = unsigned int(full_count);
				param->time_diff_s = time_diff * iFreq/full_count;

				//printf("count, time_diff, FPS >> %d %.3e %.3f\n",count,time_diff * iFreq/count,param->FPS);
				time_diff = 0.;
				local_count = 0;

				SetEvent(param->hExt_StatsUpdate);

				break;
			}

			QueryPerformanceCounter(&li_t2);
			//if (count % TP_iXon.FPSinterval == 0)
			//if ((double(li_t2.QuadPart - li_prev.QuadPart) * iFreq) >= param->display_image_interval_s)

			if (local_count >= param->display_image_count)
			{
				li_prev = li_t2;
				param->FPS = double(full_count) / (double(li_t2.QuadPart - li_t0.QuadPart) * iFreq);
		
				param->frame_count = unsigned int(full_count);
				param->time_diff_s = time_diff * iFreq/local_count;
				SetEvent(param->hExt_StatsUpdate);

				//printf("count, time_diff, FPS >> %d %.3e %.3f\n",count,time_diff * iFreq/count,param->FPS);
				time_diff = 0.;
				local_count = 0;
				localCount = 0;
			}
		}

		param->bCalibration = false;		// TURN OFF CALIBRATION FLAG

		ResetEvent(param->hStartAcquisition);
	}
	else
	{
		//while (camconf->nleft > 0)
		//while(WaitForAcquisitionTimeOut(param->acqTimeout_ms) == DRV_SUCCESS)
		while (CallWaitForAcquisitionTimeout(param) == SUCCESS)
		{
			QueryPerformanceCounter(&li_t1);
			// CODE FOR 'Internal trigger mode' MODE
			//error = GetNumberNewImages(&first_im, &last_im);
			CallGetLastNumberImages(&TP_Simu_iXon, param, first_im, last_im);

			//// CODE FOR 'Software trigger' MODE
			//SendSoftwareTrigger();
			//WaitForAcquisition();
			//GetNumberNewImages(&first_im, &last_im);

			if (param->bSimulation) {
				GetSimuLastImage(param);
			}
			else {
				error = GetMostRecentImage(&param->vImgBuffer[0], imgSize);
			}

			delta = last_im - last_read_im;
			if (delta > 1)
			{
				//cout << "delta image : " << last_im - last_read_im << last_im << last_read_im << endl;
				SetEvent(param->hExt_MissedFrame);
				SetEvent(TP_Simu_iXon.hMissedFrame);
			}
			last_read_im = last_im;
		

			if (param->bAsyncProc == true)
			{
				if (TryEnterCriticalSection(&param->CSCopyImgBuffer) == TRUE)
				{
					std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(),param->vImgBufferCopy.begin());
				}
			}
			else
			{
				// FOR TEST PURPOSE!!!!
				std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(),param->vImgBufferCopy.begin());

				// SHARED MEMORY DISPLAY IMAGE BUFFER
				//SHM_imarray[0].md[0].write = 1;				// signal you are about to write
				//std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(),(long*)SHM_imarray[0].array.SI32);

				retVal = AOC_ProcessSHData(param);
				QueryPerformanceCounter(&li_t2);

				time_diff += double(li_t2.QuadPart - li_t1.QuadPart);

				if (retVal == SUCCESS)
				{

				}
			}

			//SetEvent(param->hFrameReady);

			// CHECK FOR LOOP TERMINATION
			if (WaitForSingleObject(param->hStopAcquisition,0) == WAIT_OBJECT_0)
			{
				//camconf->nleft = 0;
				if (::AbortAcquisition() != DRV_SUCCESS)
				{
					cout << string_format("'::AbortAcquisition()' FAILED!\n");
				}
				else
				{
					cout << string_format("Acquisition stopped successfully\n\n");
				}
				param->bStreamON = false;
				break;
			}

			++full_count;
			++local_count;
			param->testLocalCount++;
			localCount++;

			QueryPerformanceCounter(&li_t2);
			//if (count % TP_iXon.FPSinterval == 0)
			//if ((double(li_t2.QuadPart - li_prev.QuadPart) * iFreq) >= param->display_image_interval_s)

			// cout << "ixon simu: " << localCount << endl;

			if (localCount >= param->display_image_count)
			{
				li_prev = li_t2;
				param->FPS = double(full_count) / (double(li_t2.QuadPart - li_t0.QuadPart) * iFreq);
				SHM_imarray[0].md[0].write = 1;				// signal you are about to write
				std::copy(param->vImgBuffer.begin(), param->vImgBuffer.end(), (long*)SHM_imarray[0].array.SI32);

				param->frame_count = unsigned int(full_count);
				param->time_diff_s = time_diff * iFreq/local_count;
				SetEvent(param->hExt_StatsUpdate);

				//printf("count, time_diff, FPS >> %d %.3e %.3f\n",count,time_diff * iFreq/count,param->FPS);
				time_diff = 0.;
				local_count = 0;
				localCount = 0;
			}
		}
	}

	//t1 = (float)(now.tv_nsec)*1e-9 + (float)(ptm->tm_sec);
	//printf("\n");
	//double avg_time = (t1 - t0)/ifr;
	//double avg_freq = 1./avg_time;
	//printf("Average time between images (freq): %.3e (%.3f)\n",avg_time,avg_freq);

	param->bAcqON		= false; // updating acquisition flags
	param->bStreamON	= false; // updating acquisition flags

	return NULL;
}

mcUINT32 Thread_iXon(void *pParam)
{
	ThreadParams_iXon *param = (ThreadParams_iXon*)pParam;

	DWORD WaitResult;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hShutdownCamera);
	vEvent.push_back(param->hStartAcquisition);
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
				// SHUTDOWN THE CAMERA
				::ShutDown();
				break;
			}
		case WAIT_OBJECT_0+2:
			{
				// START ACQUISITION
				AcquireLoop(param);
				SetEvent(param->hAcquisitionDone);
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

void printBits(unsigned int num)
{
	for(unsigned int bit=0;bit<(sizeof(unsigned int) * 8); bit++)
	{
		printf("%i ", num & 0x01);
		num = num >> 1;
	}

	printf("\n");
}

mcINT16 CameraSelect (int camID)
{
	at_32 lNumCameras;
	GetAvailableCameras(&lNumCameras);
    
	if ((camID >= 0) && (camID < lNumCameras))
	{
		at_32 lCameraHandle;
		GetCameraHandle(camID, &lCameraHandle);
		SetCurrentCamera(lCameraHandle);
		return camID;
	}
	else
	{
		return SERIOUS;
	}

	return SUCCESS;
}

mcINT16 AOC_SetHWConf(ThreadParams_iXon *param)
{
	unsigned long error;

	// ------------------------------------------
	//    startup configuration of the camera
	// ------------------------------------------
	//param->readMode = 4;
	//param->acqMode = 5;
	//param->FTmode = 0;
	//param->ampType = 0;
	//param->ADCchannel = 0;
	//param->iVHSIndex = 1;

	error = SetPCIMode(1,1);						// DO NOT CHANGE THESE VALUES (CHANGING ONLY RELEVANT WITH CCI23 PCI CARD)
	error = SetDMAParameters(1, 0.0005f);
	SetFrameTransferMode(param->FTmode);			// Frame transfer mode is OFF
	SetReadMode(param->readMode);					// Set Read Mode to --Image--
	SetAcquisitionMode(param->acqMode);				// Set Acq. mode to --Run till abort--
	int bInternalShutter = 0;
	error = IsInternalMechanicalShutter(&bInternalShutter);
	error = SetShutter(1,1,0,0);					// Initialize Shutter

	error = SetOutputAmplifier(param->ampType);
	error = SetADChannel(param->ADCchannel);

	//GetFastestRecommendedVSSpeed (&iVHSIndex,&fs);
	SetVSSpeed(param->iVHSIndex);

	error = ::SetNumberAccumulations(1);
	error = GetAcquisitionTimings(&(param->expTime_s), &(param->accTime_s), &(param->kinTime_s));
	error = GetReadOutTime(&param->readoutTime);

	TP_iXon.display_image_count = 50;//int(1. / TP_iXon.kinTime_s * TP_iXon.display_image_interval_s + 0.5);
	
	printf("For AD channel = %d: ", param->ADCchannel);
	printf("\n >> Exp. time = %.6f s, acc. time = %.6f s, kin. time = %.6f s, readout time = %.6f\n",param->expTime_s, param->accTime_s, param->kinTime_s, param->readoutTime);

	return SUCCESS;
}

mcINT16 AOC_InitializeCamera(ThreadParams_iXon* param)
{

	unsigned long error;
	int gain = -1;
	int temp;
	//int temp_min = 1;
	//int temp_max = 1;
	unsigned int state;

	int glow = -1;
	int ghigh = -1;
	int nb_pre_gains = -1;
	int nb_adcs = -1;
	int nb_amps = -1;
	int nb_speeds = -1;
	//int adc_channel = -1;
	//int amp_type = -1;
	//float pre_gain = 0.0;
	//float hzs_speed = 0.0;

	float fval = 0.f;

	param->bCamInitialized = false;

	//char *cmd_in = (char*) malloc(LINESZ * sizeof(char));
	//char *buf = (char*) malloc(LINESZ * sizeof(char));
	//char *auxstr = (char*) malloc(LINESZ * sizeof(char));
	//char *msgout = (char*) malloc(LINESZ * sizeof(char));
	//mytest *= 1;

	//cconf = (cam_config*) malloc(sizeof(cam_config));
	//cam_config_init(cconf);

	//sem_t *test_sem = (sem_t*)my_sem_open("test_semaphore",O_CREAT,0644,1);


	// ------------------------------------------
	//            Initialize CCD
	// ------------------------------------------

	// CAMERA SIMULATION MODE: -1
	if (!param->bSimulation) {
		if (CameraSelect(param->camID) < 0)
		{
			cout << "\n\nCamera working in SIMULATION mode.\n"
				"==================================\n\n";
			param->bCamInitialized = true;
			param->bSimulation = true;
			return SUCCESS;
		}
		else
		{
			cout << "\n\nCamera initialization in progress (may take a few seconds)\n"
				"==========================================================\n\n";
		}

		error = Initialize("");

		if (error != DRV_SUCCESS)
		{
			cout << "Initialisation error...exiting" << endl;
			return SERIOUS;
		}

		//// ------------------------------------------
		////    startup configuration of the camera
		//// ------------------------------------------
		param->readMode = 4;
		param->acqMode = 5;
		param->FTmode = 1;								// Frame transfer mode is OFF (0) or ON (1)

		error = SetPCIMode(1, 1);						// DO NOT CHANGE THESE VALUES (CHANGING ONLY RELEVANT WITH CCI23 PCI CARD)
		error = SetDMAParameters(1, 0.001f);
		SetFrameTransferMode(param->FTmode);			// Frame transfer mode
		SetReadMode(param->readMode);					// Set Read Mode to --Image--
		SetAcquisitionMode(param->acqMode);				// Set Acq. mode to --Run till abort--
		SetShutter(1, 1, 0, 0);							// Initialize Shutter
		SetExposureTime(param->expTime_s);				// Set initial exposure time
		GetDetector(&(param->DX), &(param->DY));		// Get Detector dimensions  
		SetImage(param->binX, param->binY, param->startX, param->DX, param->startY, param->DY);		// Setup Image dimensions

		SetKineticCycleTime(0.f); // test
		error = SetNumberAccumulations(1);


		//error = GetAcquisitionTimings(&(param->expTime_s),&(param->accTime_s), &(param->kinTime_s));
		//printf("\n >> Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s\n",param->expTime_s, param->accTime_s, param->kinTime_s);

		AndorCapabilities* caps = (AndorCapabilities*)malloc(sizeof(AndorCapabilities));
		caps->ulSize = sizeof(AndorCapabilities);
		error = GetCapabilities(caps);
		if (error != DRV_SUCCESS)
		{
			cout << "GET Capabilities fail" << endl;
			return SERIOUS;
		}

		cout << "Capabilities:" << endl << "============" << endl;
		printf("Acq Modes  : "); printBits(caps->ulAcqModes);
		printf("Read Modes : "); printBits(caps->ulReadModes);
		printf("Trig Modes : "); printBits(caps->ulTriggerModes);
		printf("Cam Type   : %u \n", caps->ulCameraType);
		printf("Pixl Mode  : "); printBits(caps->ulPixelMode);
		printf("Set Funcs  : "); printBits(caps->ulSetFunctions);
		printf("Get Funcs  : "); printBits(caps->ulGetFunctions);
		printf("Get Feats  : "); printBits(caps->ulFeatures);
		printf("PCI speed  : %u \n", caps->ulPCICard);
		printf("EMG capab  : "); printBits(caps->ulEMGainCapability);
		printf("FT Readmod : "); printBits(caps->ulFTReadModes);
		cout << "============" << endl;

		////ensure frame transfer is off
		//SetFrameTransferMode (0);
		SetBaselineClamp(0);

		error = GetNumberPreAmpGains(&nb_pre_gains);
		printf("# of Preamp Gains : %d\n", nb_pre_gains);

		for (int i = 0; i < nb_pre_gains; i++)
		{
			error = GetPreAmpGain(i, &fval);
			printf("%d: gain = %.2f\n", i, fval);
		}
		param->PreAmpGainIndex = 0;
		error = SetPreAmpGain(param->PreAmpGainIndex); // set pre-amp gain to 1.0

		// ====================================================================
		error = GetNumberAmp(&nb_amps);
		printf("%d output amplifier available \n", nb_amps);

		// ====================================================================
		error = GetNumberADChannels(&nb_adcs);
		printf("\n%d A/D channels are available \n", nb_adcs);

		int depth = -1;
		for (int i = 0; i < nb_adcs; i++)
		{
			error = GetBitDepth(i, &depth);
			printf("- channel #%d: depth = %d\n", i, depth);
		}
		printf("\n\n");
		// ====================================================================

		//amp_type    = 0; // type of output amplification (1: conventional, 0: EMCCD)
		//adc_channel = 0; //
		param->ampType = 0;
		param->ADCchannel = 0;

		error = SetOutputAmplifier(param->ampType);
		error = SetADChannel(param->ADCchannel);

		//printf("For AD channel = %d: ", param->ADCchannel);
		error = GetNumberHSSpeeds(param->ADCchannel, param->ampType, &nb_speeds);
		printf("%d speeds are available in this mode\n", nb_speeds);

		for (int i = 0; i < nb_speeds; i++)
		{
			error = GetHSSpeed(param->ADCchannel, param->ampType, i, &fval);
			printf("- i=%d: speed = %.1f MHz\n", i, fval);
		}
		printf("\n");

		// ====================================================================

		//SetHighCapacity(1); // 
		error = GetEMCCDGain(&gain);

		if (error != DRV_SUCCESS)
			cout << "GET EMCCDGain fail" << endl;

		printf("EMCCD Gain = %d\n", gain);

		GetEMGainRange(&glow, &ghigh);
		printf("EM Gain range = %d - %d\n", glow, ghigh);

		GetMCPGain(&gain);
		printf("MCP   Gain = %d\n", gain);

		GetMCPGainRange(&glow, &ghigh);
		printf("MCP Gain range = %d - %d\n", glow, ghigh);

		param->expTime_s = 0.00001f;
		//param->kinTime_s = 0.0f;
		//SetExposureTime(param->expTime_s);
		//SetKineticCycleTime(param->kinTime_s);
		////float readoutTime = -1.f;
		//param->readoutTime = -1.f;

		//error = GetReadOutTime(&param->readoutTime);
		//error = GetAcquisitionTimings(&(param->expTime_s),&(param->accTime_s), &(param->kinTime_s));
		//
		//printf("\n >> Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s, readout time = %.6f\n",param->expTime_s, param->accTime_s, param->kinTime_s, param->readoutTime);

		AOC_SetHWConf(&TP_iXon);

		param->trigMode = 0;				// 0: Internal trigger | 10: Software trigger
		SetTriggerMode(param->trigMode);
		if (IsTriggerModeAvailable(param->trigMode) != DRV_SUCCESS)
		{
			printf("Trigger mode %d is not available!\n", param->trigMode);
		}

		if (caps->ulTriggerModes & AC_TRIGGERMODE_CONTINUOUS)
		{
			printf("Correct Hardware for Continuous Mode.\n");
		}
		else
		{
			printf("InCorrect Hardware for Continuous Mode.\n");
		}

		param->iVHSIndex = 1;
		//GetFastestRecommendedVSSpeed (&iVHSIndex,&fs);
		SetVSSpeed(param->iVHSIndex);

		//float VSspeed = -1.f;
		GetVSSpeed(param->iVHSIndex, &param->VSspeed);
		printf("Vertical Speed set to %g microseconds per pixel shift\r\n", param->VSspeed);

		if (caps->ulSetFunctions & AC_SETFUNCTION_BASELINECLAMP)
		{
			error = SetBaselineClamp(1);
			if (error != DRV_SUCCESS)
			{
				printf("Set Baseline Clamp Error\n");
			}
		}

		param->vImgBufferCopy.resize(param->DX * param->DY);

		// ------------------------------------------
		//    special setup for the temperature
		// ------------------------------------------
		if (false)
		{
			state = CoolerON();
			cout << "CoolerON status code: " << state << endl;

			state = GetTemperatureRange(&param->tempMin, &param->tempMax);
			cout << "Detector temperature range = " << param->tempMin << " - " << param->tempMax << endl;

			state = GetTemperature(&temp);    // get detector temperature
			cout << "Detector temperature = " << temp << " deg" << endl;

			temp = 0; // temperature setting (air cooling)
			state = SetTemperature(temp);

			state = GetTemperature(&temp);    // get detector temperature
			cout << "Detector temperature = " << temp << " deg" << endl;

			cout << "GetTemperature status code: " << state << endl;

			while (state != DRV_TEMPERATURE_STABILIZED)
			{
				Sleep(1000);
				state = GetTemperature(&temp);    // get detector temperature
				printf("\rTemperature = %+03d deg - status code: %6d", temp, state);
				fflush(stdout);
			}

			switch (state)
			{
			case DRV_TEMPERATURE_OFF:
				cout << "Cooler OFF" << endl; break;
			case DRV_TEMPERATURE_STABILIZED:
				cout << "Stabilized" << endl; break;
			case DRV_TEMPERATURE_NOT_REACHED:
				cout << "Cooling" << endl; break;
			default:
				cout << "Cooler status unknown" << endl;
			}
		}
		else
		{
			float fTemp = -999.f;
			GetTemperatureF(&fTemp);
			state = ::CoolerOFF();
			param->bCoolingON = false;
			cout << "CoolerOFF status code: " << state << endl;
			cout << "Current temperature: " << fTemp << endl;
		}

			param->bCamInitialized = true;
	}
	else
	{
		param->expTime_s = 0.001;
		param->DX = TP_calcSH.DX;
		param->DY = TP_calcSH.DY;
		param->vImgBufferCopy.resize(param->DX* param->DY);
		param->bCamInitialized = true;
	}
	return SUCCESS;
}

mcINT16 AOC_ShutdownCamera(ThreadParams_iXon *param)
{
	if (param->bCamInitialized == false)
	{
		return WARNING;
	}

	if (param->bSimulation == true)
	{
		param->bCamInitialized = false;
		return SUCCESS;
	}

	int error = ::ShutDown();

	if (error != DRV_SUCCESS)
	{
		LogMsg = string_format("[%s] FAILED to shutdown the camera.",__FUNCTION__);
		cout << LogMsg << endl;
		Log->StackAddMsg(LOG_SERIOUS,0,__LINE__,__FILEREF__,LogMsg);
		return SERIOUS;
	}

	return SUCCESS;
}

mcINT16 GetCS_CopyImgBuffer(CRITICAL_SECTION *_CS)
{
	_CS = &(TP_iXon.CSCopyImgBuffer);

	return SUCCESS;
}

mcINT16 AOC_ProcessSHData(ThreadParams_iXon *param)
{
	//// COMPUTE MIN AND MAX IMAGE VALUES FOR FAST MEDIAN COMPUTATION
	//long min = 0;
	//long max = 0;
	//param->sci.MinMax(param->vImgBufferCopy,min,max);
	//// COMPUTE FAST MEDIAN
	//double bkg = param->sci.Median(param->vImgBufferCopy,min,max);

	double bkg = 0.;

	// USED FOR THRESHOLD DETERMINATION ONLY ==> THRESHOLD SETTING CAN BE DONE BY PYTHON CODE (SLOW/ARBITRARY/USER CHANGES)
	if (true)
	{
		// COMPUTE MEAN
		vector<double> mom;
		param->sci.Moment(param->vImgBufferCopy, 1, mom);
		bkg = mom[0];
	}
	else
	{
		bkg = 0.;
	}

	// COMPUTE MULTI-THREADED IMAGE CENTROIDS
	trigger_calc_SH_data_mt(&param->vImgBufferCopy[0], float(bkg));

	return SUCCESS;
}