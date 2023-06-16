#include <cstdlib>
#include <vector>
#include <iostream>

#include <MCODE.h>

#include "ImageStruct.h"
#include "thread_DataLogger_params.h"
#include "thread_calcSH_params.h"
#include "thread_DMcmd_params.h"

#undef SUCCESS
#include <asdkDM.h>

using namespace std;
using namespace MultiCODE;

extern MultiCODE::MCODE_LogManager	*Log;

extern int					N_DM_COMMAND_MAPS;

extern ThreadParams_DataLogger	TP_DataLogger;
extern DWORD 					Thread_DataLogger_ID;
extern HANDLE					hThread_DataLogger;

extern ThreadParams_calcSH	TP_calcSH;
extern DWORD 				Thread_calcSH_ID;
extern HANDLE				hThread_calcSH;

static int					nThreads = 0;
static HANDLE				hThreadManager = INVALID_HANDLE_VALUE;
static int					ThreadManagerID = 0;

extern acs::DM				dm;
extern vector<acs::Scalar>	vDM_cmd;

extern vector<float>		vtmpLO;
extern vector<float>		vtmpHO;
extern vector<double>		vTT_drift;
extern vector<double>		vHO_drift;

extern float				leakage_gain;

extern vector<float>		vCmdMatrix;

extern IMAGE				*SHM_imarray;			// SHARED MEMORY INSTANT IXON IMAGE
extern IMAGE				*SHM_Slopes;			// SHARED MEMORY N FRAMES HISTORY DATA (DISPLAY RATE)
extern IMAGE				*SHM_RefCent;			// SHARED MEMORY REFERENCE CENTROIDS DATA
extern IMAGE				*SHM_DM_slope_mask;		// SHARED MEMORY PUPIL MASK DATA
extern IMAGE				*SHM_DM_state;			// SHARED MEMORY DM STATE
extern IMAGE				*SHM_ModesInfo;			// SHARED MEMORY MODE INFORMATION
extern IMAGE				*SHM_Timestamps;        // SHARED MEMORY TIMESTAMPS
extern vector<IMAGE*>		vSHM_DM_BUFFERS;		// SHARED MEMORY

extern vector<vector<float>>	vDM_cmd_buffers;		// DM CORRECTION BUFFERS (N_DM_COMMAND_MAPS x nAct)
extern vector<float>			vDM_static_cmd;			// DM STATIC COMMAND MAP (MAP CORRESPONDING TO NON-TIME-CRITICAL MAPS)
extern vector<size_t>			vCollapseStaticInd;		// DM STATIC MAPS INDICES
extern vector<size_t>			vCollapseDynamicInd;	// DM DYNAMIC (TIME-CRITICAL) INDICES = NOT(vCollapseStaticInd)
extern vector<float>			vDM_cur_cmd;			// DM INSTANTANEOUS COMMAND MAP = vDM_cmd_buffers[vCollapseIndices] + vDM_static_cmd
extern vector<long long>		vTimeStamp;				// TIME VECTOR RECORDED JUST BEFORE SETTING THE DM INTO POSITION

extern vector<int>				vDMSquare2DMmap;

	// VALID SUBAPERTURE INDEX
extern vector<int>				vSH_ValidCell;

extern int						localCount;

MCODE_Export_FITS FITS;


LARGE_INTEGER t0;
LARGE_INTEGER t1;

LARGE_INTEGER liFreq;


static inline int calc_SH_data_mt(int threadID, int sizex, int sizey, long *vImg, float mean_bkg, int nx, int ny, int nsub, int *v_indx, int *v_indy, 
	int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp, float *SH_slope_x, float *SH_slope_y, float *SH_phot, float *SH_refx, float *SH_refy, float *v_thresh, int *v_enabled)
{
	//cout << "sizex: " << sizex << endl;
	//cout << "sizey: " << sizey << endl;
	//cout << "nx   : " << nx << endl;
	//cout << "ny   : " << ny << endl;

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
	//double thresh = double(1.3 * mean_bkg);
	//double thresh = double(mean_bkg) + 100.;

	//cout << "[" << threadID << "]" << " thresh: " << thresh << endl;

	int xs;
	int ys;

	int ii;
	int jj;

	//LARGE_INTEGER t0;
	//LARGE_INTEGER t1;

	//LARGE_INTEGER liFreq;
	//QueryPerformanceFrequency(&liFreq);
	//QueryPerformanceCounter(&t0);

	//printf("[%d] nsub: %d %x %x %x %x %x %x %x\n", threadID, nsub, v_xx, v_yy, v_indx, v_indy, SH_xtmp, SH_ytmp, SH_phot);

	int write = 0;

	for (int i = 0; i < nsub; ++i)
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


		// FIXME: PRE-COMPUTE SUBAPERTURE INDICES
		// EXTRACT THE SUBPUPIL DATA FROM ORIGINAL IMAGE
		sub_arr.resize(xs * ys);

		int ycrd;
		int ycrd_sa;
		//#pragma omp parallel for
		for (int v = 0; v < size_t(ys); ++v)
		{
			ycrd = (y0 + v) * sizex;
			ycrd_sa = v * xs;
			for (int u = 0; u < size_t(xs); ++u)
			{
				sub_arr[u + ycrd_sa] = vImg[x0 + u + ycrd];
			}
		}

		/*std::string imagettes = "imagettes" + std::to_string(i);
		FITS.Save(&sub_arr[0], xs, ys, "R:\\", imagettes);*/


		/*for (int a = 0; a != 100; a++) {
			std::string imagettes = "imagettes" + to_string(a);
			FITS.Save(&sub_arr[a], , 13, "R:\\", imagettes);
			a++;
		}*/

		// FIND THE MAXIMUM VALUE IN THE SUBPUPIL
		maxval = -INT_MAX;
		vector<int>::iterator ptr = sub_arr.begin();

		//#pragma omp parallel for shared (maxval)
		for (int k = 0; k < xs * ys; ++k)
		{
			if (*ptr > maxval) maxval = *ptr;
			++ptr;
		}

		// LINEAR INDEX OF THE CURRENT SH APERTURE
		int crd = ii + jj * nx;

		float xc = 0.f;
		float yc = 0.f;
		float sum = 0.f;

		//printf("[%d] maxval >> %d %.3f\n",threadID,maxval,thresh);
		//SH_phot[crd] = float(maxval);

		//================================
		// FOR DEBUGGING PURPOSE ONLY!!
		//thresh = 0.;
		//================================

		//if (maxval > thresh)
		if ((v_enabled[i] == 1) && (maxval > v_thresh[i]))
		{
			// SUBTRACT MIN
//#pragma omp parallel for
			for (int k = 0; k < xs * ys; ++k)
			{
				//if (sub_arr[k] < thresh) sub_arr[k] = 0;
				if (sub_arr[k] < v_thresh[i]) sub_arr[k] = 0;
				// OR
				// if (sub_arr[k] < thresh) sub_arr[k] -= thresh;
			}

			// COMPUTE CENTROID
			xc = 0.f;
			yc = 0.f;
			sum = 0.;
			float val;

			//#pragma omp parallel for reduction(+:xc, yc, sum)
			for (register int v = 0; v < ys; ++v)
			{
				ycrd_sa = v * xs;
				for (register int u = 0; u < xs; ++u)
				{
					val = float(sub_arr[u + ycrd_sa]);
					xc += val * float(u);
					yc += val * float(v);
					sum += val;
				}
			}

			//cout << "sum >> " << sum << endl;

			if (sum > 0.f)
			{
				xc = xc / sum;
				yc = yc / sum;
				//cout << "xc, yc >> " << xc << "  " << yc << endl;
			}
			else
			{
				xc = float(0.5f * (xs - 1));
				yc = float(0.5f * (ys - 1));
			}

		}
		else
		{
			// SET TO REF POSITION
			xc = float(0.5f * (xs - 1));
			yc = float(0.5f * (ys - 1));
			sum = 0.f;
		}

		//printf("xc,yc >> %.3f %.3f\n",xc,yc);

		//SH_xtmp[crd] = float(xc) - SH_refx[crd];
		//SH_ytmp[crd] = float(yc) - SH_refy[crd];
		SH_xtmp[crd] = float(xc);
		SH_ytmp[crd] = float(yc);
		SH_slope_x[crd] = SH_xtmp[crd] - SH_refx[crd];
		SH_slope_y[crd] = SH_ytmp[crd] - SH_refy[crd];
		//cout << "crd, sh_slope_x sh_slope_y: " << crd << " " << SH_slope_x[crd] << " " << SH_slope_y[crd] << endl;
		SH_phot[crd] = float(sum);
		//cout << "image " << i << " x y " << xc << " " << yc << " " << SH_slope_x[crd] << " " << SH_slope_y[crd] << endl;
		/*if (write == 0) {
			ofstream file;
			file.open("C:\\Users\\lucas\\Documents\\STAGE\\Misc\\data_aoc_slopes.txt", fstream::app);
			if (!file) {
				cout << "cant open file" << endl;
			}
			else {
				file << SH_slope_x[crd] << " " << SH_slope_y[crd] << endl;
			}
			file.close();

		}
		if (i >= 99)
			write = 1;*/
	}
	//QueryPerformanceCounter(&t1);
	//printf("Execution time [%d]: %.3e\n",threadID, double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart));

	return 0;
}

static DWORD Thread_CalcSHData(void *pParam)
{
	ThreadParams_CalcSHData *param = (ThreadParams_CalcSHData*)pParam;


	vector<HANDLE> vEvent;
	vEvent.push_back(param->hTerminateThread);
	vEvent.push_back(param->hStartProc);
	int nEvents = int(vEvent.size());

	SetEvent(param->hThreadCreated);

	//printf("Thread %d, param->counter = %d\n",param->threadNumber,param->counter);
	//printf("Thread %d, param->counter = %d  param: %x\n",param->threadNumber,param->counter, param);

	//LARGE_INTEGER t0;
	//LARGE_INTEGER t1;

	//LARGE_INTEGER liFreq;
	//QueryPerformanceFrequency(&liFreq);

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
				//QueryPerformanceCounter(&t0);
				calc_SH_data_mt(
					param->threadNumber,
					param->DX,
					param->DY,
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
					param->SH_slp_x,
					param->SH_slp_y,
					param->SH_phot,
					param->SH_refx,
					param->SH_refy,
					param->vThresh,
					param->vEnabled);
				param->counter++;

				//QueryPerformanceCounter(&t1);
				//printf("Execution time [%d]: %.3e\n",param->threadNumber, double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart));

				SetEvent(param->hProcDone);
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

	return 0;
}

mcINT16 set_thread_params(int threadID, int _nsub, vector<int>& _v_xx, vector<int>& _v_yy, vector<int>& _vInd_i, vector<int>& _vInd_j, 
	vector<float>& _SH_xtmp, vector<float>& _SH_ytmp, vector<float>& _SH_slope_x, vector<float>& _SH_slope_y, vector<float>& _SH_phot, 
	vector<float>& _SH_refx, vector<float>& _SH_refy, vector<float>& _vThresh, vector<int>& _vEnabled)
{
	if ((threadID < 0) || (threadID >= TP_calcSH.nThreads))
	{
		cout << string_format("[%s] invalid thread ID '%d' (0...%d).\n",__FUNCTION__,threadID,TP_calcSH.nThreads-1);
		return WARNING;
	}

	if (TP_calcSH.logVerbose > 0)
	{
		cout << string_format("v_xx/yy index/y [%d]: %x %x %x %x  %d %d\n",threadID,&_v_xx[0],&_v_yy[0],&_vInd_i[0],&_vInd_j[0], _vInd_i[0],_vInd_j[0]);
		cout << string_format("SH_xtmp/y [%d]: %x %x %x\n",threadID,&_SH_xtmp[0],&_SH_ytmp[0],&_SH_phot[0]);
		cout << string_format("vInd_i: [");
		for (size_t k=0;k<_nsub;++k)
		{
			cout << string_format(" %d",_vInd_i[k]);
		}
		cout << string_format("]\nvInd_j: [");
		for (size_t k=0;k<_nsub;++k)
		{
			cout << string_format(" %d",_vInd_j[k]);
		}
		cout << string_format("]\n");
	}

	TP_calcSH.vTP_CalcSHData[threadID].nsub = _nsub;
	TP_calcSH.vTP_CalcSHData[threadID].v_xx = &_v_xx[0];
	TP_calcSH.vTP_CalcSHData[threadID].v_yy = &_v_yy[0];
	TP_calcSH.vTP_CalcSHData[threadID].vIndx = &_vInd_i[0];
	TP_calcSH.vTP_CalcSHData[threadID].vIndy = &_vInd_j[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_xtmp = &_SH_xtmp[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_ytmp = &_SH_ytmp[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_slp_x = &_SH_slope_x[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_slp_y = &_SH_slope_y[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_phot = &_SH_phot[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_refx = &_SH_refx[0];
	TP_calcSH.vTP_CalcSHData[threadID].SH_refy = &_SH_refy[0];
	size_t index = size_t(TP_calcSH.nx*TP_calcSH.ny / TP_calcSH.nThreads * threadID);
	TP_calcSH.vTP_CalcSHData[threadID].vThresh = &_vThresh[index];
	TP_calcSH.vTP_CalcSHData[threadID].vEnabled = &_vEnabled[index];
	return 0;
}

mcINT16 SetSubApertureThresholds(int _nThreads, vector<float>& _vThresh)
{

	int nCells = TP_calcSH.nx * TP_calcSH.ny;
	int nsub = nCells / _nThreads;

	size_t count = 0;

	for (int k = 0; k < _nThreads - 1; ++k)
	{
		for (size_t u=0;u<nsub;++u)
		{
			TP_calcSH.vThresh[count] = _vThresh[count];
			++count;
		}
	}

	int rem_sub = nCells - nsub * (TP_calcSH.nThreads - 1);
	for (size_t u = 0; u < rem_sub; ++u)
	{
		TP_calcSH.vThresh[count] = _vThresh[count];
		++count;
	}

	return 0;
}

mcINT16 StartThread_calcSH(int _nThreads, int _DX, int _DY, int _nx, int _ny)
{
	TP_calcSH.hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_calcSH.hTerminateThread = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_calcSH.hStartProc = CreateEvent(NULL,FALSE,FALSE,NULL);
	TP_calcSH.hProcDone = CreateEvent(NULL,FALSE,FALSE,NULL);

	TP_calcSH.hOut_SingleModeDone = CreateEvent(NULL,FALSE,FALSE,"_AOC_CALIB_SINGLE_MODE_DONE");
	TP_calcSH.hOut_FullCalibDone = CreateEvent(NULL,FALSE,FALSE,"_AOC_CALIB_FULL_CALIB_DONE");

	TP_calcSH.nThreads = _nThreads;
	nThreads = _nThreads;

	printf("Trying to create %d threads!\n",TP_calcSH.nThreads);

	TP_calcSH.rad_per_pix_x = 1.;
	TP_calcSH.rad_per_pix_y = 1.;

	TP_calcSH.bEnableCloseLoop = false;

	TP_calcSH.bCalDataAvailable_ZER = false;
	TP_calcSH.nModes_ZER = 10;
	TP_calcSH.vDM_Modes_ZER_Slope.resize(TP_calcSH.nModes_ZER * TP_calcSH.nx*TP_calcSH.ny * 2);
	TP_calcSH.vDM_Modes_ZER_iNorm.resize(TP_calcSH.nModes_ZER);
	TP_calcSH.vDM_ProjOnModes_ZER.resize(TP_calcSH.nModes_ZER);

	TP_calcSH.bCalDataAvailable_ZON = false;
	TP_calcSH.nModes_ZON = 97;
	TP_calcSH.vDM_Modes_ZON_Slope.resize(TP_calcSH.nModes_ZON * TP_calcSH.nx*TP_calcSH.ny * 2);
	TP_calcSH.vDM_Modes_ZON_iNorm.resize(TP_calcSH.nModes_ZON);
	TP_calcSH.vDM_ProjOnModes_ZON.resize(TP_calcSH.nModes_ZON);

	TP_calcSH.vTP_CalcSHData.resize(TP_calcSH.nThreads);
	TP_calcSH.vThreadHandles.resize(TP_calcSH.nThreads);
	TP_calcSH.vThresh.resize(TP_calcSH.nThreads);
	TP_calcSH.vThread_CalcSHDataID.resize(TP_calcSH.nThreads);

	TP_calcSH.timestamp.resize(50);

	vector<HANDLE> vCreationEvents(TP_calcSH.nThreads);

	DWORD waitResult;

	for (size_t k=0;k<TP_calcSH.nThreads;++k)
	{
		TP_calcSH.vTP_CalcSHData[k].bEndThread = false;
		TP_calcSH.vTP_CalcSHData[k].hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
		TP_calcSH.vTP_CalcSHData[k].hTerminateThread = CreateEvent(NULL,FALSE,FALSE,NULL);
		TP_calcSH.vTP_CalcSHData[k].hStartProc = CreateEvent(NULL,FALSE,FALSE,NULL);
		TP_calcSH.vTP_CalcSHData[k].hProcDone = CreateEvent(NULL,FALSE,FALSE,NULL);

		vCreationEvents[k] = TP_calcSH.vTP_CalcSHData[k].hThreadCreated;

		TP_calcSH.vTP_CalcSHData[k].nx = _nx;
		TP_calcSH.vTP_CalcSHData[k].ny = _ny;
		TP_calcSH.vTP_CalcSHData[k].nsub = 0;
		TP_calcSH.vTP_CalcSHData[k].SH_xtmp = NULL;
		TP_calcSH.vTP_CalcSHData[k].SH_ytmp = NULL;
		TP_calcSH.vTP_CalcSHData[k].SH_slp_x = NULL;
		TP_calcSH.vTP_CalcSHData[k].SH_slp_y = NULL;
		TP_calcSH.vTP_CalcSHData[k].SH_phot = NULL;
		TP_calcSH.vTP_CalcSHData[k].SH_refx = NULL;
		TP_calcSH.vTP_CalcSHData[k].SH_refy = NULL;
		TP_calcSH.vTP_CalcSHData[k].DX = _DX;
		TP_calcSH.vTP_CalcSHData[k].DY = _DY;
		TP_calcSH.vTP_CalcSHData[k].vImg = NULL;
		TP_calcSH.vTP_CalcSHData[k].bkg = 0.f;
		TP_calcSH.vTP_CalcSHData[k].vIndx = NULL;
		TP_calcSH.vTP_CalcSHData[k].vIndy = NULL;
		TP_calcSH.vTP_CalcSHData[k].v_xx = NULL;
		TP_calcSH.vTP_CalcSHData[k].v_yy = NULL;

		TP_calcSH.vTP_CalcSHData[k].counter = 0;
		TP_calcSH.vTP_CalcSHData[k].threadNumber = int(k);

		TP_calcSH.vThreadHandles[k] = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_CalcSHData,&TP_calcSH.vTP_CalcSHData[k],THREAD_PRIORITY_HIGHEST,&TP_calcSH.vThread_CalcSHDataID[k]);

		cout << string_format("HANDLE for thread %d: %x  %ud\n",k,TP_calcSH.vThreadHandles[k],TP_calcSH.vThread_CalcSHDataID[k]);
	}

	if (TP_calcSH.nThreads > 0)
	{
		waitResult = WaitForMultipleObjects(DWORD(TP_calcSH.vTP_CalcSHData.size()),&vCreationEvents[0],TRUE,THREAD_CREATION_TIMEOUT);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			cout << string_format("All %d threads created successfully!\n",TP_calcSH.nThreads);
			break;
		case WAIT_TIMEOUT:
			cout << string_format("The creation of %d threads FAILED! (TIMEOUT)\n",TP_calcSH.nThreads);
			break;
		}
	}

	return 0;
}

mcINT16 TerminateThread_calcSH()
{
	DWORD waitResult;

	if (TP_calcSH.vTP_CalcSHData.size() != 0)
	{
		vector<HANDLE> vEvents(TP_calcSH.vTP_CalcSHData.size());

		for (size_t k=0;k<TP_calcSH.vTP_CalcSHData.size();++k)
		{
			SetEvent(TP_calcSH.vTP_CalcSHData[k].hTerminateThread);
		}

		waitResult = WaitForMultipleObjects(DWORD(TP_calcSH.vTP_CalcSHData.size()),&TP_calcSH.vThreadHandles[0],TRUE,THREAD_TERMINATION_TIMEOUT);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			printf("All %d threads successfully terminated!\n",TP_calcSH.nThreads);
			break;
		case WAIT_TIMEOUT:
			printf("TIEMOUT while trying to terminate all %d threads.\n",TP_calcSH.nThreads);
			break;
		}
	}

	return 0;
}

int trigger_calc_SH_data_mt(long *vImg, float _bkg)
{
	int nEvents;
	vector<HANDLE> vWaitEvent;
	LARGE_INTEGER now;

	QueryPerformanceFrequency(&liFreq);

	for (size_t k=0;k<nThreads;++k)
	{
		vWaitEvent.push_back(TP_calcSH.vTP_CalcSHData[k].hProcDone);
		TP_calcSH.vTP_CalcSHData[k].vImg = vImg;
		TP_calcSH.vTP_CalcSHData[k].bkg = _bkg;
	}

	nEvents = int(vWaitEvent.size());

	// TRIGGER THE SLOPES COMPUTATION
	for (size_t k=0;k<nThreads;++k)
	{
		SetEvent(TP_calcSH.vTP_CalcSHData[k].hStartProc);
	}

	bool bTriggerDisplay = false;

	LARGE_INTEGER li_t0;
	LARGE_INTEGER li_t1;
	double pc_freq = 0.;
	double pc_start = 0.;

	double tick_period_microsecond = 1. / double(liFreq.QuadPart) * 1.e6;

	QueryPerformanceCounter(&li_t0);

	int waitResult = WaitForMultipleObjects(nEvents,&vWaitEvent[0],TRUE,5000);
	// AT THIS POINT, ALL THE SLOPES HAVE BEEN FRESHLY COMPUTED

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
	{
		int nSubAp = TP_calcSH.nx*TP_calcSH.ny;
		int mod_counter = localCount % TP_calcSH.display_image_count;
		QueryPerformanceCounter(&li_t1);

		double t = ((li_t1.QuadPart - li_t0.QuadPart) * tick_period_microsecond) / liFreq.QuadPart;

		TP_calcSH.timestamp[mod_counter] = li_t1.QuadPart;

		if (mod_counter != 0) {
			double test_allez = (TP_calcSH.timestamp[mod_counter] - TP_calcSH.timestamp[mod_counter - 1]) / 10000000;
			//cout << "valeur test de la diff: " << test_allez << endl;
		}


		//cout << timestamp[mod_counter] << " " << mod_counter << endl;
		// ===============================================================================
		// ===============================================================================
		// ============================================== PUT THE SLOPE ANALYSIS CODE HERE

		//cout << "calcSH: " << localCount << "and mod counter: " << mod_counter << endl;

		if (TP_calcSH.bCalMode == false)
		{
			// COMPUTE PROJECTION OF CURRENT MAP ONTO THE CALIBRATED MODES
			switch (TP_calcSH.projMode)
			{
			case TP_calcSH.AOC_PROJ_MODE_ZER:
				trigger_proj_ZER();
				std::memcpy(&TP_calcSH.vModesDataHist[TP_calcSH.MaxModes * mod_counter], &TP_calcSH.vDM_ProjOnModes_ZER[0], sizeof(float) * TP_calcSH.vDM_ProjOnModes_ZER.size());

				// GET THE TMESTAMP AND SET THE DM INTO POSITION!!

				QueryPerformanceCounter(&now);
				trigger_set_DM_correction_map();
				vTimeStamp[mod_counter] = now.QuadPart;
				break;
			case TP_calcSH.AOC_PROJ_MODE_ZON:
				trigger_proj_ZON();
				std::memcpy(&TP_calcSH.vModesDataHist[TP_calcSH.MaxModes * mod_counter], &TP_calcSH.vDM_ProjOnModes_ZON[0], sizeof(float) * TP_calcSH.vDM_ProjOnModes_ZON.size());
				break;
			}
		}

		// ===============================================================================
		// ===============================================================================
		// ===============================================================================

		// COPY EACH SUBAPERTURE INSTANT SLOPE DATA (vSH_cent_x, vSH_cent_y, vSH_slope_x, vSH_slope_y, vSH_phot) INTO THE FLAT vFlatSHData VECTOR
		//std::memcpy(&TP_calcSH.vFlatSHData[TP_calcSH.nx*TP_calcSH.ny*mod_counter*TP_calcSH.nDataSets]		,&TP_calcSH.vSH_cent_x[0],sizeof(float)*nSubAp);
		//std::memcpy(&TP_calcSH.vFlatSHData[TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 1)]	,&TP_calcSH.vSH_cent_y[0],sizeof(float)*nSubAp);
		//std::memcpy(&TP_calcSH.vFlatSHData[TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 2)]	,&TP_calcSH.vSH_slope_x[0],sizeof(float)*nSubAp);
		//std::memcpy(&TP_calcSH.vFlatSHData[TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 3)]	,&TP_calcSH.vSH_slope_y[0],sizeof(float)*nSubAp);
		//std::memcpy(&TP_calcSH.vFlatSHData[TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 4)]	,&TP_calcSH.vSH_phot[0],sizeof(float)*nSubAp);

		size_t bufNum = TP_calcSH.circBufCounter % TP_calcSH.nCircBuffers;
		std::memcpy(&TP_calcSH.vCircBuf_FlatSHData[bufNum][TP_calcSH.nx*TP_calcSH.ny*mod_counter*TP_calcSH.nDataSets]		,&TP_calcSH.vSH_cent_x[0],sizeof(float)*nSubAp);
		std::memcpy(&TP_calcSH.vCircBuf_FlatSHData[bufNum][TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 1)]	,&TP_calcSH.vSH_cent_y[0],sizeof(float)*nSubAp);
		std::memcpy(&TP_calcSH.vCircBuf_FlatSHData[bufNum][TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 2)]	,&TP_calcSH.vSH_slope_x[0],sizeof(float)*nSubAp);
		std::memcpy(&TP_calcSH.vCircBuf_FlatSHData[bufNum][TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 3)]	,&TP_calcSH.vSH_slope_y[0],sizeof(float)*nSubAp);
		std::memcpy(&TP_calcSH.vCircBuf_FlatSHData[bufNum][TP_calcSH.nx*TP_calcSH.ny*(mod_counter*TP_calcSH.nDataSets + 4)]	,&TP_calcSH.vSH_phot[0],sizeof(float)*nSubAp);
		std::memcpy(&TP_calcSH.vCircBuf_Timestamps[bufNum][50 * mod_counter], &TP_calcSH.timestamp[0], sizeof(double)*50);


		TP_DataLogger.circBufferCalcSHCounter = bufNum;
		SetEvent(TP_DataLogger.hDataBlockAvailable);

		//cout << TP_calcSH.local_counter % TP_calcSH.display_image_count << " " << mod_counter << endl;

		++TP_calcSH.local_counter;
		
		
		if (TP_calcSH.bCalMode == false)
		{
			bTriggerDisplay = mod_counter != 0 ? false : true;
		}
		else
		{
			//TP_calcSH.local_counter = (++TP_calcSH.local_counter) % TP_calcSH.nCalFrameCount;

			vector<float>::iterator it;

			// FILL THE CORRESPONDING VECTOR ACCORDING TO THE TYPE OF CALIBRATION (ZER, ZON)
			if (TP_calcSH.bCalDataAvailable_ZER == true)
			{
				it = TP_calcSH.vDM_Modes_ZER_Slope.begin() + TP_calcSH.curCalMode*nSubAp*2;
			}
			else if (TP_calcSH.bCalDataAvailable_ZON == true)
			{
				it = TP_calcSH.vDM_Modes_ZON_Slope.begin() + TP_calcSH.curCalMode*nSubAp*2;
			}

			// ACCUMULATE THE X/Y SLOPE VECTORS
			for (size_t k=0;k<nSubAp;++k,++it)
			{
				(*it) += TP_calcSH.vSH_slope_x[k];
			}
			for (size_t k=0;k<nSubAp;++k,++it)
			{
				(*it) += TP_calcSH.vSH_slope_y[k];
			}

			bTriggerDisplay = mod_counter != 0 ? false : true;

			//if ((localCount % TP_calcSH.nCalFrameCount) == 0)
			//{
			//	SetEvent(TP_calcSH.hOut_SingleModeDone);
			//}
		}

		if (bTriggerDisplay == true)
		{
			// WHEN THE NUMBER OF ACQUIRED DATA SETS EQUALS TP_calcSH.display_image_c	ount THEN TRANSFER THE DATA BLOCK TO THE SHARED MEMORY 
			//std::copy(TP_calcSH.vFlatSHData.begin(),TP_calcSH.vFlatSHData.end(),(float*)SHM_Slopes[0].array.F);
			
			//timestamp.clear();
			//std::copy(TP_calcSH.timestamp.begin(), TP_calcSH.timestamp.end(), (long long *)SHM_Timestamps[0].array.D);
			std::copy(TP_calcSH.vCircBuf_Timestamps[bufNum].begin(), TP_calcSH.vCircBuf_Timestamps[bufNum].end(), SHM_Timestamps[0].array.D);
			std::copy(TP_calcSH.vCircBuf_FlatSHData[bufNum].begin(),TP_calcSH.vCircBuf_FlatSHData[bufNum].end(),(float*)SHM_Slopes[0].array.F);
			std::copy(TP_calcSH.vModesDataHist.begin(),TP_calcSH.vModesDataHist.end(),(float*)SHM_ModesInfo[0].array.F);
			std::copy(vDM_cur_cmd.begin(),vDM_cur_cmd.end(),(float*)SHM_DM_state[0].array.F);
			++TP_calcSH.circBufCounter;
		}
		break;
	}
	case WAIT_TIMEOUT:
		printf("TIMEOUT!!\n");
		break;
	}

	return 0;
}

int trigger_proj_TT()
{


	return 0;
}

int trigger_proj_ZER()
{
	if (false)
	{
		vector<float>::iterator itRecModes_x;
		vector<float>::iterator itRecModes_y;
		vector<float>::iterator itSlopes_x;
		vector<float>::iterator itSlopes_y;

		// PROJECT MEASURED SLOPES ONTO THE SLOPES OF THE CALIBRATED MODES
		for (size_t i = 0; i < TP_calcSH.nModes_ZER; ++i)
		{
			itRecModes_x = TP_calcSH.vDM_Modes_ZER_Slope.begin() + 2 * i * TP_calcSH.nx * TP_calcSH.ny;
			itRecModes_y = TP_calcSH.vDM_Modes_ZER_Slope.begin() + (2 * i + 1) * TP_calcSH.nx * TP_calcSH.ny;
			itSlopes_x = TP_calcSH.vSH_slope_x.begin();
			itSlopes_y = TP_calcSH.vSH_slope_y.begin();

			float val = 0.f;

			for (size_t k = 0; k < TP_calcSH.nx * TP_calcSH.ny; ++k, ++itRecModes_x, ++itSlopes_x, ++itRecModes_y, ++itSlopes_y)
			{
				val += (*itSlopes_x) * (*itRecModes_x) + (*itSlopes_y) * (*itRecModes_y);
			}

			if (TP_calcSH.vDM_Modes_ZER_iNorm[i] != 0.f)
			{
				val *= TP_calcSH.vDM_Modes_ZER_iNorm[i];
			}

			// THIS VALUE CONTAINS THE AMPLITUDE OF EACH OF THE CALIBRATED MODES CONTAINED IN THE MEASURED SLOPE MAP
			TP_calcSH.vDM_ProjOnModes_ZER[i] = val;
		}
	}

	if (true)
	{
		vector<float>::iterator itCmdMat_x;
		vector<float>::iterator itCmdMat_y;
		vector<float>::iterator itSlopes_x;
		vector<float>::iterator itSlopes_y;

		// USE THE COMMAND MATRIX TO COMPUTE THE COMMAND COEFFICIENTS
		for (size_t i = 0; i < TP_calcSH.nModes_ZER; ++i)
		{
			itCmdMat_x = vCmdMatrix.begin() + 2 * i * TP_calcSH.nx * TP_calcSH.ny;
			itCmdMat_y = vCmdMatrix.begin() + (2 * i + 1) * TP_calcSH.nx * TP_calcSH.ny;
			itSlopes_x = TP_calcSH.vSH_slope_x.begin();
			itSlopes_y = TP_calcSH.vSH_slope_y.begin();

			float val = 0.f;

			for (size_t k = 0; k < TP_calcSH.nx * TP_calcSH.ny; ++k, ++itCmdMat_x, ++itSlopes_x, ++itCmdMat_y, ++itSlopes_y)
			{
				val += (*itSlopes_x) * (*itCmdMat_x) + (*itSlopes_y) * (*itCmdMat_y);
			}

			//if (TP_calcSH.vDM_Modes_ZER_iNorm[i] != 0.f)
			//{
			//	val *= TP_calcSH.vDM_Modes_ZER_iNorm[i];
			//}

			// THIS VALUE CONTAINS THE AMPLITUDE OF EACH OF THE CALIBRATED MODES CONTAINED IN THE MEASURED SLOPE MAP
			TP_calcSH.vDM_ProjOnModes_ZER[i] = val;
		}
	}

	return 0;
}

int trigger_set_DM_correction_map()
{
	if (TP_calcSH.bEnableCloseLoop == true)
	{
		// TIP TILT
		for (size_t i = 0; i < 2; ++i)
		{
			vTT_drift[i] = vTT_drift[i]*leakage_gain + TP_calcSH.ZER_gain * (-TP_calcSH.vDM_ProjOnModes_ZER[i]);

			for (size_t k = 0; k < TP_calcSH.nAct; ++k)
			{
				vtmpLO[k] = vtmpLO[k]*leakage_gain + TP_calcSH.ZER_gain * (-TP_calcSH.vDM_ProjOnModes_ZER[i] * TP_calcSH.vDM_Cmd_ZER[i * TP_calcSH.nAct + k]);
			}
		}

		vSHM_DM_BUFFERS[2]->kw[0].value.numf = vTT_drift[0];
		vSHM_DM_BUFFERS[2]->kw[1].value.numf = vTT_drift[1];

		std::copy(vtmpLO.begin(), vtmpLO.end(), vSHM_DM_BUFFERS[2]->array.F);
		vSHM_DM_BUFFERS[2]->md[0].cnt0++;

		if (true)
		{
			// HIGHER ORDER
			for (size_t i = 2; i < TP_calcSH.nModes_ZER; ++i)
			{
				vHO_drift[i-2] = vHO_drift[i-2]*leakage_gain + TP_calcSH.ZER_gain * (-TP_calcSH.vDM_ProjOnModes_ZER[i]);

				for (size_t k = 0; k < TP_calcSH.nAct; ++k)
				{
					vtmpHO[k] = vtmpHO[k]*leakage_gain + TP_calcSH.ZER_gain * (-TP_calcSH.vDM_ProjOnModes_ZER[i] * TP_calcSH.vDM_Cmd_ZER[i * TP_calcSH.nAct + k]);
				}
			}

			for (size_t i = 0; i < TP_calcSH.nModes_ZER - 2; ++i)
			{
				vSHM_DM_BUFFERS[3]->kw[i].value.numf = vHO_drift[i];
			}

			std::copy(vtmpHO.begin(), vtmpHO.end(), vSHM_DM_BUFFERS[3]->array.F);
			vSHM_DM_BUFFERS[3]->md[0].cnt0++;
		}

		collapse_DM_buffers();
		DM_set_to_position();
	}
	return 0;
}

int trigger_proj_ZON()
{

	return 0;
}

