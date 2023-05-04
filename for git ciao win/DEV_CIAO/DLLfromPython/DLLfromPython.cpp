#include <Windows.h>

#include <iostream>
#include <vector>

#include "DLLfromPython.h"

using namespace std;

#define THREAD_TERMINATION_TIMEOUT	5000
#define	THREAD_CREATION_TIMEOUT		5000


typedef struct
{
	bool		bEndThread;
	HANDLE		hThreadCreated;
	HANDLE		hTerminatedThread;
	HANDLE		hStartProc;
	HANDLE		hProcDone;

	int threadNumber;
	int counter;

	int			sizex;
	int			sizey;
	int*		vImg;
	float		bkg;
	int			nx;
	int			ny;
	int			nsub;
	int*		vIndx;
	int*		vIndy;
	int*		v_xx;
	int*		v_yy;
	float*		SH_xtmp;
	float*		SH_ytmp;
	float*		SH_phot;
} ThreadParams_CalcSHData;


typedef struct
{
	bool			bEndThread;
	HANDLE			hThreadCreated;
	HANDLE			hTerminatedThread;
	HANDLE			hInitialize;
	HANDLE			hStartMonitor;
	HANDLE			hStopMonitor;
	HANDLE			hProcDone;

	bool			bMonitorInitialized;

	// EXTERNAL EVENTS
	int				nEvents;
	vector<HANDLE>	vTimingEvent;
	vector<int>		vEventIndex;

	vector<double>	vTimingSum;
	int				counter;

	int				displayFreq;
} ThreadParams_MonitorTimings;

static int			nThreads = 0;
static HANDLE		hThreadManager = INVALID_HANDLE_VALUE;
static int			ThreadManagerID = 0;

//static int			counter = 0;

static bool			bEndThread = false;
static HANDLE		hThreadCreated = INVALID_HANDLE_VALUE;
static HANDLE		hTerminatedThread = INVALID_HANDLE_VALUE;
static HANDLE		hStartProc = INVALID_HANDLE_VALUE;
static HANDLE		hProcDone = INVALID_HANDLE_VALUE;

vector<HANDLE>		vThreadHandles;
vector<ThreadParams_CalcSHData>
					vTP_CalcSHData;
vector<DWORD>		vThread_CalcSHDataID;


ThreadParams_MonitorTimings	TP_MonitorTimings;
DWORD Thread_MonitorTimingsID;
HANDLE hThread_MonitorTimings;

extern "C"
int nop()
//int calc_SH_data(int argc, void *argv[])
{
	return 0;
}

extern "C"
int calc_SH_data(int sizex, int sizey, int *vImg, float mean_bkg, int nx, int ny, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp, float *SH_phot)
{
	int ncx = nx;
	int ncy = ny;

	//cout << "sizex: " << sizex << endl;
	//cout << "sizey: " << sizey << endl;
	//cout << "nx   : " << nx << endl;
	//cout << "ny   : " << ny << endl;

	//int bkg_ui = 0;
	//for (int k=0;k<sizex*sizey;++k)
	//{
	//	bkg_ui += vImg[k];
	//}
	//float mean_bkg = float(bkg_ui)/float(sizex*sizey);

	int x0;
	int y0;
	int x1;
	int y1;

	vector<int> sub_arr;

	int maxval;
	double thresh = double(1.3 * mean_bkg);

	//cout << "thresh: " << thresh << endl;

	int xs;
	int ys;

	//LARGE_INTEGER t0;
	//LARGE_INTEGER t1;

	//LARGE_INTEGER liFreq;
	//QueryPerformanceFrequency(&liFreq);
	//QueryPerformanceCounter(&t0);

	for (int j=0;j<ncy;++j)
	{
		y0 = v_yy[j];
		y1 = v_yy[j + 1];

		for (int i=0;i<ncx;++i)
		{
			x0 = v_xx[i];
			x1 = v_xx[i + 1];

			xs = x1-x0;
			ys = y1-y0;

			//cout << "x0,x1,y0,y1,xs,ys: " << x0 << "," << x1 << "," << y0 << "," << y1 << "," << xs << "," << ys << endl;

			// EXTRACT THE SUBPUPIL DATA FROM ORIGINAL IMAGE
			sub_arr.resize(xs*ys);

			for (size_t v=0;v<size_t(ys);++v)
			{
				for (size_t u=0;u<size_t(xs);++u)
				{
					sub_arr[u + v*xs] = vImg[x0+u + (y0+v)*sizex];
				}
			}
			
			// FIND THE MAXIMUM VALUE IN THE SUBPUPIL
			maxval = -INT_MAX;
			vector<int>::iterator ptr = sub_arr.begin();
			for (int k=0;k<xs*ys;++k, ++ptr)
			{
				if (*ptr > maxval) maxval = *ptr;
			}

			int crd = i + j*nx;

			float xc = 0.f;
			float yc = 0.f;
			float sum = 0.f;

			SH_phot[crd] = float(maxval);

			if (maxval > thresh)
			{
				// SUBTRACT MIN
				for (int k=0;k<xs*ys;++k)
				{
					if (sub_arr[k] < thresh) sub_arr[k] = 0;
					// OR
					// if (sub_arr[k] < thresh) sub_arr[k] -= thresh;
				}

				// COMPUTE CENTROID
				xc = 0.f;
				yc = 0.f;
				sum = 0.;
				float val;

				for (register int i=0;i<xs;++i)
				{
					for (register int j=0;j<ys;++j)
					{
						val = float(sub_arr[i + j*xs]);
						xc += val * float(i);
						yc += val * float(j);
						sum += val;
					}
				}

				if (sum > 0.f)
				{
					xc = xc / sum;
					yc = yc / sum;
				}
				else
				{
					xc = float(0.5f * (xs-1));
					yc = float(0.5f * (ys-1));
				}

			}
			else
			{
				// SET TO REF POSITION
				xc = float(0.5f * (xs-1));
				yc = float(0.5f * (ys-1));
			}

			//SH_xtmp[crd] = float(x0 + xc);
			//SH_ytmp[crd] = float(y0 + yc);
			SH_xtmp[crd] = float(xc);
			SH_ytmp[crd] = float(yc);
		}
	}

	//QueryPerformanceCounter(&t1);

	//cout << "Execution time: " << double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart) << endl;

	return 0;
}

// sizex		[IN][int]		x size of the input vImg array (pix)
// sizey		[IN][int]		y size of the input vImg array (pix)
// *vImg		[IN][int*]		pointer to the image array
// nx			[IN][int]		total number of subapertures of the WFS in x
// ny			[IN][int]		total number of subapertures of the WFS in y
// n_sub		[IN][int]		number of subapertures x coordinates to be processed bu this thread (also the size of the v_indx, and v_indy arrays)
// *v_indx		[IN][int]		x index array of the subapertures
// *v_indy		[IN][int]		y index array of the subapertures
// *v_xx		[IN][int]		pointer to the i-th subaperture pixel 'x' coordinate
// *v_yy		[IN][int]		pointer to the i-th subaperture pixel 'y' coordinate
// *SH_xtmp		[OUT][float]	pointer to the photometric barycenter x coordinate of the i-th subaperture
// *SH_ytmp		[OUT][float]	pointer to the photometric barycenter y coordinate of the i-th subaperture
// *SH_phot		[OUT][float]	pointer to the photometric illumination of the i-th subaperture
//int calc_SH_data_mt(int sizex, int sizey, int *vImg, float mean_bkg, int nx, int ny, int nsub, int *v_indx, int *v_indy, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp, float *SH_phot)
int calc_SH_data_mt(int threadID, int sizex, int sizey, int *vImg, float mean_bkg, int nx, int ny, int nsub, int *v_indx, int *v_indy, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp, float *SH_phot)
{
	//cout << "sizex: " << sizex << endl;
	//cout << "sizey: " << sizey << endl;
	//cout << "nx   : " << nx << endl;
	//cout << "ny   : " << ny << endl;

	//// FIXME:: COMPUTE THE MEAN BACKGROUND OUTSIDE EACH THREAD!
	//int bkg_ui = 0;
	//for (int k=0;k<sizex*sizey;++k)
	//{
	//	bkg_ui += vImg[k];
	//}
	//float mean_bkg = float(bkg_ui)/float(sizex*sizey);
	//printf("bkg, vImg, nsub >> %.3f %d %d\n",mean_bkg, vImg[5 + 5*10], nsub);

	//for (size_t k=0;k<nsub;++k)
	//{
	//	printf("xx,yy: %d %d\n",v_indx[k],v_indy[k]);
	//}

	int x0;
	int y0;
	int x1;
	int y1;

	vector<int> sub_arr;

	int maxval;
	double thresh = double(1.3 * mean_bkg);

	//cout << "thresh: " << thresh << endl;

	int xs;
	int ys;

	//LARGE_INTEGER t0;
	//LARGE_INTEGER t1;

	//LARGE_INTEGER liFreq;
	//QueryPerformanceFrequency(&liFreq);
	//QueryPerformanceCounter(&t0);

	int ii;
	int jj;

	//printf("[%d] nsub: %d %x %x %x %x %x %x %x\n", threadID, nsub, v_xx, v_yy, v_indx, v_indy, SH_xtmp, SH_ytmp, SH_phot);

	for (int i=0;i<nsub;++i)
	{
		ii = v_indx[i];
		jj = v_indy[i];

		//printf("[%d] ii,jj >> %d %d\n",threadID,ii,jj);

		y0 = v_yy[jj];
		y1 = v_yy[jj + 1];

		x0 = v_xx[ii];
		x1 = v_xx[ii + 1];

		//printf("x0,x1, y0,y1 >> %d %d -- %d %d\n",x0,x1,y0,y1);

		xs = x1 - x0;
		ys = y1 - y0;

		//cout << "x0,x1,y0,y1,xs,ys: " << x0 << "," << x1 << "," << y0 << "," << y1 << "," << xs << "," << ys << endl;

		// EXTRACT THE SUBPUPIL DATA FROM ORIGINAL IMAGE
		sub_arr.resize(xs*ys);

		for (size_t v=0;v<size_t(ys);++v)
		{
			for (size_t u=0;u<size_t(xs);++u)
			{
				sub_arr[u + v*xs] = vImg[x0+u + (y0+v)*sizex];
			}
		}
			
		// FIND THE MAXIMUM VALUE IN THE SUBPUPIL
		maxval = -INT_MAX;
		vector<int>::iterator ptr = sub_arr.begin();
		for (int k=0;k<xs*ys;++k, ++ptr)
		{
			if (*ptr > maxval) maxval = *ptr;
		}

		int crd = ii + jj*nx;

		float xc = 0.f;
		float yc = 0.f;
		float sum = 0.f;

		//printf("[%d] maxval >> %d %.3f\n",threadID,maxval,thresh);
		SH_phot[crd] = float(maxval);

		if (maxval > thresh)
		{
			// SUBTRACT MIN
			for (int k=0;k<xs*ys;++k)
			{
				if (sub_arr[k] < thresh) sub_arr[k] = 0;
				// OR
				// if (sub_arr[k] < thresh) sub_arr[k] -= thresh;
			}

			// COMPUTE CENTROID
			xc = 0.f;
			yc = 0.f;
			sum = 0.;
			float val;

			for (register int i=0;i<xs;++i)
			{
				for (register int j=0;j<ys;++j)
				{
					val = float(sub_arr[i + j*xs]);
					xc += val * float(i);
					yc += val * float(j);
					sum += val;
				}
			}

			if (sum > 0.f)
			{
				xc = xc / sum;
				yc = yc / sum;
			}
			else
			{
				xc = float(0.5f * (xs-1));
				yc = float(0.5f * (ys-1));
			}

		}
		else
		{
			// SET TO REF POSITION
			xc = float(0.5f * (xs-1));
			yc = float(0.5f * (ys-1));
		}

		//printf("xc,yc >> %.3f %.3f\n",xc,yc);

		SH_xtmp[crd] = float(xc);
		SH_ytmp[crd] = float(yc);
	}

	//QueryPerformanceCounter(&t1);

	//printf("Execution time [%d]: %.3e\n",threadID, double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart));

	return 0;
}

static DWORD Thread_CalcSHData(void *pParam)
{
	ThreadParams_CalcSHData *param = (ThreadParams_CalcSHData*)pParam;


	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminatedThread);
	vEvent.push_back(param->hStartProc);
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	//printf("Thread %d, param->counter = %d\n",param->threadNumber,param->counter);
	//printf("Thread %d, param->counter = %d  param: %x\n",param->threadNumber,param->counter, param);

	int waitResult;

	while(param->bEndThread == false)
	{
		waitResult = WaitForMultipleObjects(nEvents,&vEvent[0],FALSE,INFINITE);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			{
				param->bEndThread = true;
				break;
			}
		case WAIT_OBJECT_0+1:
			{
				//printf("Current thread [%d] STARTED! %x\n",param->threadNumber,param->vImg);
				calc_SH_data_mt(
					param->threadNumber,
					param->sizex,
					param->sizey,
					param->vImg,
					param->bkg,
					param->nx,
					param->ny,
					param->nsub,
					param->vIndx,
					param->vIndy,
					param->v_xx,
					param->v_yy,
					param->SH_xtmp,
					param->SH_ytmp,
					param->SH_phot);
				param->counter++;
				//printf("Current thread (%d) counter value: %d\n",param->threadNumber,param->counter);
				SetEvent(param->hProcDone);
				break;
			}
		case WAIT_FAILED:
			{
				printf("WAIT_FAILED!!\n");
				break;
			}
		case WAIT_ABANDONED:
			{
				printf("WAIT_ABANDONED!!\n");
				break;
			}
		}
	}

	return 0;
}

extern "C"
int StartThread_ThreadManager(int _nThreads, int _sizex, int _sizey, int _nx, int _ny)
{
	hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	hTerminatedThread = CreateEvent(NULL,FALSE,FALSE,NULL);
	hStartProc = CreateEvent(NULL,FALSE,FALSE,NULL);
	hProcDone = CreateEvent(NULL,FALSE,FALSE,NULL);

	nThreads = _nThreads;

	printf("Trying to create %d threads!\n",nThreads);

	DWORD waitResult;

	vTP_CalcSHData.resize(nThreads);
	vThreadHandles.resize(nThreads);
	vThread_CalcSHDataID.resize(nThreads);
	vector<HANDLE> vCreationEvents(nThreads);

	for (size_t k=0;k<nThreads;++k)
	{
		vTP_CalcSHData[k].bEndThread = false;
		vTP_CalcSHData[k].hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
		vTP_CalcSHData[k].hTerminatedThread = CreateEvent(NULL,FALSE,FALSE,NULL);
		vTP_CalcSHData[k].hStartProc = CreateEvent(NULL,FALSE,FALSE,NULL);
		vTP_CalcSHData[k].hProcDone = CreateEvent(NULL,FALSE,FALSE,NULL);

		vCreationEvents[k] = vTP_CalcSHData[k].hThreadCreated;

		vTP_CalcSHData[k].nx = _nx;
		vTP_CalcSHData[k].ny = _ny;
		vTP_CalcSHData[k].nsub = 0;
		vTP_CalcSHData[k].SH_xtmp = NULL;
		vTP_CalcSHData[k].SH_ytmp = NULL;
		vTP_CalcSHData[k].SH_phot = NULL;
		vTP_CalcSHData[k].sizex = _sizex;
		vTP_CalcSHData[k].sizey = _sizey;
		vTP_CalcSHData[k].vImg = NULL;
		vTP_CalcSHData[k].bkg = 0.f;
		vTP_CalcSHData[k].vIndx = NULL;
		vTP_CalcSHData[k].vIndy = NULL;
		vTP_CalcSHData[k].v_xx = NULL;
		vTP_CalcSHData[k].v_yy = NULL;

		vTP_CalcSHData[k].counter = 0;
		vTP_CalcSHData[k].threadNumber = int(k);

		vThreadHandles[k] = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_CalcSHData,&vTP_CalcSHData[k],THREAD_PRIORITY_NORMAL,&vThread_CalcSHDataID[k]);

		printf("HANDLE for thread %d: %x  %ud\n",k,vThreadHandles[k],vThread_CalcSHDataID[k]);
	}

	if (nThreads > 0)
	{
		waitResult = WaitForMultipleObjects(DWORD(vTP_CalcSHData.size()),&vCreationEvents[0],TRUE,THREAD_CREATION_TIMEOUT);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			printf("All %d threads created successfully!\n",nThreads);
			break;
		case WAIT_TIMEOUT:
			printf("The creation of %d threads FAILED!\n",nThreads);
			break;
		}
	}

	return 0;
}

extern "C"
int TerminateThread_ThreadManager(int nThreads)
{
	DWORD waitResult;

	if (vTP_CalcSHData.size() != 0)
	{
		vector<HANDLE> vEvents(vTP_CalcSHData.size());

		for (size_t k=0;k<vTP_CalcSHData.size();++k)
		{
			SetEvent(vTP_CalcSHData[k].hTerminatedThread);
		}

		waitResult = WaitForMultipleObjects(DWORD(vTP_CalcSHData.size()),&vThreadHandles[0],TRUE,THREAD_TERMINATION_TIMEOUT);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			printf("All %d threads successfully terminated!\n",nThreads);
			break;
		case WAIT_TIMEOUT:
			printf("TIEMOUT while trying to terminate all %d threads.\n",nThreads);
			break;
		}
	}

	return 0;
}

extern "C"
int trigger_calc_SH_data(int sizex, int sizey, int *vImg, float bkg, int nx, int ny, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp, float *SH_phot)
{
	return calc_SH_data(sizex, sizey, vImg, bkg, nx, ny, v_xx, v_yy, SH_xtmp, SH_ytmp, SH_phot);
}

extern "C"
int set_thread_params(int threadID, int _nsub, int *_v_xx, int *_v_yy, int *_vInd_i, int *_vInd_j, float *_SH_xtmp, float *_SH_ytmp, float *_SH_phot)
{
	if ((threadID < 0) || (threadID >= nThreads))
	{
		return -1;
	}

	printf("v_xx/yy index/y [%d]: %x %x %x %x  %d %d\n",threadID,_v_xx,_v_yy,_vInd_i,_vInd_j, _vInd_i[0],_vInd_j[0]);
	printf("SH_xtmp/y [%d]: %x %x %x\n",threadID,_SH_xtmp,_SH_ytmp,_SH_phot);
	printf("vInd_i DLL: [");
	for (size_t k=0;k<_nsub;++k)
	{
		printf(" %d",_vInd_i[k]);
	}
	printf("]\nvInd_j DLL: [");
	for (size_t k=0;k<_nsub;++k)
	{
		printf(" %d",_vInd_j[k]);
	}
	printf("]\n");

	vTP_CalcSHData[threadID].nsub = _nsub;
	vTP_CalcSHData[threadID].v_xx = _v_xx;
	vTP_CalcSHData[threadID].v_yy = _v_yy;
	vTP_CalcSHData[threadID].vIndx = _vInd_i;
	vTP_CalcSHData[threadID].vIndy = _vInd_j;
	vTP_CalcSHData[threadID].SH_xtmp = _SH_xtmp;
	vTP_CalcSHData[threadID].SH_ytmp = _SH_ytmp;
	vTP_CalcSHData[threadID].SH_phot = _SH_phot;

	return 0;
}

extern "C"
//int trigger_calc_SH_data_mt(int sizex, int sizey, int *vImg, int nx, int ny, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp)
int trigger_calc_SH_data_mt_test()
{
	int nEvents;
	vector<HANDLE> vWaitEvent;

	for (size_t k=0;k<nThreads;++k)
	{
		vWaitEvent.push_back(vTP_CalcSHData[k].hProcDone);
	}

	nEvents = int(vWaitEvent.size());

	for (int iter=0;iter<2;++iter)
	{
		for (size_t k=0;k<nThreads;++k)
		{
			SetEvent(vTP_CalcSHData[k].hStartProc);
		}

		int waitResult = WaitForMultipleObjects(nEvents,&vWaitEvent[0],TRUE,5000);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			printf("All done!!\n");
			break;
		case WAIT_TIMEOUT:
			printf("TIMEOUT!!\n");
			break;
		}

		Sleep(1000);
	}

	return 0;
}

extern "C"
//int trigger_calc_SH_data_mt(int sizex, int sizey, int *vImg, int nx, int ny, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp)
int trigger_calc_SH_data_mt(int *vImg, float _bkg)
{
	int nEvents;
	vector<HANDLE> vWaitEvent;

	for (size_t k=0;k<nThreads;++k)
	{
		vWaitEvent.push_back(vTP_CalcSHData[k].hProcDone);
		vTP_CalcSHData[k].vImg = vImg;
		vTP_CalcSHData[k].bkg = _bkg;
	}

	nEvents = int(vWaitEvent.size());

	for (size_t k=0;k<nThreads;++k)
	{
		SetEvent(vTP_CalcSHData[k].hStartProc);
	}

	int waitResult = WaitForMultipleObjects(nEvents,&vWaitEvent[0],TRUE,5000);

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
		//printf("All done!!\n");
		break;
	case WAIT_TIMEOUT:
		printf("TIMEOUT!!\n");
		break;
	}

	return 0;
}





static DWORD Thread_MonitorTimings(void *pParam)
{
	ThreadParams_MonitorTimings *param = (ThreadParams_MonitorTimings*)pParam;

	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminatedThread);
	vEvent.push_back(param->hInitialize);
	vEvent.push_back(param->hStartMonitor);
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	//printf("Thread %d, param->counter = %d\n",param->threadNumber,param->counter);
	//printf("Thread %d, param->counter = %d  param: %x\n",param->threadNumber,param->counter, param);

	int waitResult;

	LARGE_INTEGER liFreq;
	QueryPerformanceFrequency(&liFreq);
	double fact = 1. / double(liFreq.QuadPart);

	int nEventsToMonitor = 0;

	vector<LARGE_INTEGER> vTimeStamps;

	unsigned long counter = 0;

	SetEvent(param->hThreadCreated);

	while(param->bEndThread == false)
	{
		waitResult = WaitForMultipleObjects(nEvents,&vEvent[0],FALSE,INFINITE);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			{
				param->bEndThread = true;
				break;
			}
		case WAIT_OBJECT_0+1:
			{
				// hInitialize

				nEventsToMonitor = int(param->vTimingEvent.size());
				cout << "nEventsToMonitor: " << nEventsToMonitor << endl;

				if (nEventsToMonitor >= 2)
				{
					vTimeStamps.resize(nEventsToMonitor);
					param->vTimingSum.resize(nEventsToMonitor);
					std::fill(param->vTimingSum.begin(),param->vTimingSum.end(),0.);
					param->counter = 0;

					param->bMonitorInitialized = true;
					

					//for(size_t k=0;k<nEventsToMonitor;++k)
					//{
					//	WaitForSingleObject(param->vTimingEvent[param->vEventIndex[k]],INFINITE);
					//	QueryPerformanceCounter(&vTimeStamps[k]);
					//}

					//for(size_t k=0;k<nEventsToMonitor-1;++k)
					//{
					//	param->vTimingSum[k] += double(vTimeStamps[k+1].QuadPart) - double(vTimeStamps[k].QuadPart);
					//	param->vCounts[k]++;
					//}

					//cout << "Overall timing >> " << (double(vTimeStamps[nEventsToMonitor-1].QuadPart) - double(vTimeStamps[0].QuadPart))*fact << endl;
				}
				else
				{
					vTimeStamps.clear();
					param->bMonitorInitialized = false;
				}

				SetEvent(param->hProcDone);

				break;
			}
		case WAIT_OBJECT_0+2:
			{
				// hStartMonitoring

				cout << "Starting to monitor timing events...\n";
				SetEvent(param->hProcDone);

				param->counter = 0;

				// RESET ALL CHRONOLICALLY EXPECTED FIRED EVENTS (IF NOT COULD LEAD TO DEATHTRAP)
				for(size_t k=1;k<nEventsToMonitor;++k)
				{
					param->vTimingEvent[k];
				}

				while (true)
				{
					if (WaitForSingleObject(param->hStopMonitor,0) == WAIT_OBJECT_0)
					{
						break;
					}

					// WAIT SUCCESSIVELY FOR EVENTS...
					for(size_t k=0;k<nEventsToMonitor;++k)
					{
						WaitForSingleObject(param->vTimingEvent[k],INFINITE);
						QueryPerformanceCounter(&vTimeStamps[k]);

						cout << "Event " << k << " fired! PC: " << double(vTimeStamps[k].QuadPart) << endl;
					}

					// COMPUTE CONSECUTIVE TIMINGS. INDEX '0' TIMING IS THE OVERALL TIMING
					param->vTimingSum[0] += (double(vTimeStamps[nEventsToMonitor-1].QuadPart) - double(vTimeStamps[0].QuadPart)) * fact;
					for(size_t k=1;k<nEventsToMonitor;++k)
					{
						param->vTimingSum[k] += (double(vTimeStamps[k].QuadPart) - double(vTimeStamps[k-1].QuadPart)) * fact;
					}
					param->counter++;

					cout << "param->counter " << param->counter << "  param->displayFreq: " << param->displayFreq << endl;
					if ((param->counter % param->displayFreq) == 0)
					{
						for(size_t k=1;k<nEventsToMonitor;++k)
						{
							cout << "Average consecutive timings [" << k+1 << "] >> " << param->vTimingSum[k] / double(param->counter) << endl;
						}
						cout << "Average overall timing >> " << param->vTimingSum[0] / double(param->counter) << endl;

						param->counter = 0;
						std::fill(param->vTimingSum.begin(),param->vTimingSum.end(),0.);
					}
				}

				SetEvent(param->hProcDone);
				break;
			}
		case WAIT_FAILED:
			{
				printf("WAIT_FAILED!!\n");
				break;
			}
		case WAIT_ABANDONED:
			{
				printf("WAIT_ABANDONED!!\n");
				break;
			}
		}
	}

	return 0;
}

extern "C"
int StartThread_MonitorTimings()
{
	TP_MonitorTimings.bEndThread = false;
	TP_MonitorTimings.hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_MonitorTimings.hTerminatedThread = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_MonitorTimings.hInitialize = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_MonitorTimings.hStartMonitor = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_MonitorTimings.hStopMonitor = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_MonitorTimings.hProcDone = CreateEvent(NULL,FALSE,FALSE,NULL);

	DWORD waitResult;

	hThread_MonitorTimings = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_MonitorTimings,&TP_MonitorTimings,THREAD_PRIORITY_NORMAL,&Thread_MonitorTimingsID);

	waitResult = WaitForMultipleObjects(DWORD(vTP_CalcSHData.size()),&TP_MonitorTimings.hThreadCreated,TRUE,THREAD_CREATION_TIMEOUT);

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
		printf("Thread_MonitorTimings thread created successfully!\n",nThreads);
		break;
	case WAIT_TIMEOUT:
		printf("The creation of Thread_MonitorTimings thread FAILED!\n",nThreads);
		break;
	}

	return 0;
}

extern "C"
int TerminateThread_MonitorTimings()
{
	DWORD waitResult;

	SetEvent(TP_MonitorTimings.hTerminatedThread);

	waitResult = WaitForSingleObject(hThread_MonitorTimings,THREAD_TERMINATION_TIMEOUT);

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
		printf("Thread_MonitorTimings thread successfully terminated!\n",nThreads);
		break;
	case WAIT_TIMEOUT:
		printf("TIEMOUT while trying to terminate Thread_MonitorTimings thread.\n",nThreads);
		break;
	}

	return 0;
}

extern "C"
int SetMonitorFrequency(int num)
{
	TP_MonitorTimings.displayFreq = num;

	return 0;
}

extern "C"
int ClearMonitoringEvents()
{
	TP_MonitorTimings.vTimingEvent.clear();

	return 0;
}

extern "C"
int AddMonitoringEvent()
{
	HANDLE newEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_MonitorTimings.vTimingEvent.push_back(newEvent);

	cout << "Added timing event." << endl;

	return 0;
}

extern "C"
int InitializeTimingEvents()
{
	SetEvent(TP_MonitorTimings.hInitialize);

	int waitResult = WaitForSingleObject(TP_MonitorTimings.hProcDone, 5000);

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
		{
			if (TP_MonitorTimings.bMonitorInitialized == false)
			{
				return -1;
			}

			cout << "InitializeTimingEvents success." << endl;

			break;
		}
	case WAIT_TIMEOUT:
		{
			cout << "ERROR: wait TIMEOUT!" << endl;

			return -1;
			break;
		}
	}

	return 0;
}

extern "C"
int StartMonitor(int flag)
{
	cout << "flag = " << flag << endl;
	if (flag == 1)
	{
		SetEvent(TP_MonitorTimings.hStartMonitor);
	}
	else if (flag == 0)
	{
		SetEvent(TP_MonitorTimings.hStopMonitor);
	}

	int waitResult = WaitForSingleObject(TP_MonitorTimings.hProcDone, 5000);

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
		{
			if (TP_MonitorTimings.bMonitorInitialized == false)
			{
				return -1;
			}

			cout << "StartMonitor success." << endl;

			break;
		}
	case WAIT_TIMEOUT:
		{
			cout << "ERROR: StartMonitor wait TIMEOUT!" << endl;

			return -1;
			break;
		}
	}

	return 0;
}


extern "C"
int TriggerTimingEvent(size_t EventNumber)
{
	if (EventNumber > TP_MonitorTimings.nEvents - 1)
	{
		cout << "ERROR! invalid event index." << endl;
		return -1;
	}

	SetEvent(TP_MonitorTimings.vTimingEvent[EventNumber]);
	return 0;
}

