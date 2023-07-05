#pragma once

#include <cstdlib>

//#include <Windows.h>
#include <cstdlib>
#include <Windows.h>

#include <MCODE.h>

using namespace MultiCODE;

typedef struct
{
	bool		bEndThread;
	HANDLE		hThreadCreated;
	HANDLE		hTerminateThread;
	HANDLE		hStartProc;
	HANDLE		hProcDone;

	int			threadNumber;
	int			counter;

	int			DX;
	int			DY;
	long*		vImg;
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
	float*		SH_slp_x;
	float*		SH_slp_y;
	float*		SH_phot;
	float*		SH_refx;
	float*		SH_refy;
	float*		vThresh;
	int*		vEnabled;
} ThreadParams_CalcSHData;


typedef struct
{
	// TRIGGER EVENTS
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;
	HANDLE			hStartProc;
	HANDLE			hProcDone;

	// EXTERNAL EVENTS TRIGGERS/SIGNALING
	HANDLE			hOut_SingleModeDone;
	HANDLE			hOut_FullCalibDone;

	bool			bEndThread;

	int				logVerbose;

	// INTERNALS
	int				nThreads;		// NUMBER OF THREADS DEDICATED TO THE SH CALCULUS
	int				sizex;			// NUMBER OF X PIXELS FROM THE SH SENSOR
	int				sizey;			// NUMBER OF Y PIXELS FROM THE SH SENSOR
	int				nx;				// NUMBER OF SH SUBAPERTURES IN 'X' (10)
	int				ny;				// NUMBER OF SH SUBAPERTURES IN 'Y' (10)
	int				DX;				
	int				DY;
	int				nact_lin_x;		// LINEAR NUMBER OF ACTUATORS IN 'X' (11)
	int				nact_lin_y;		// LINEAR NUMBER OF ACTUATORS IN 'Y' (11)
	int				nAct;			// TOTAL NUMBER OF ACTUATORS (97)
	int				nDM_buffers;	// NUMBER OF DM MAPS TO COLLAPSE BEFORE SENDING TO THE DM
	int				SH_x0;
	int				SH_y0;
	double			fSH_dx;
	double			fSH_dy;

	double			rad_per_pix_x;	// SH PIXEL SCALE IN RADIANS PER PIXEL SHIFT IN 'X' DIRECTION 
	double			rad_per_pix_y;	// SH PIXEL SCALE IN RADIANS PER PIXEL SHIFT IN 'Y' DIRECTION 

	// CELL BOUNDARIES
	vector<int>		vSH_x;			// VECTOR CONTAINING THE 'X' BOUNDARIES OF THE SH CELLS
	vector<int>		vSH_y;			// VECTOR CONTAINING THE 'Y' BOUNDARIES OF THE SH CELLS
	
	vector<int>		vSH_ValidCell;		// 10*10 VECTOR SPECIFYING WHETHER A SUBAPERTURE IS VALID OR NOT (OR SOMETHING ELSE)

	// CELL INDICES (FOR THREADS) AND CELL INDEX MAPPING
	vector<vector<int>> vInd_i;			// [MT] nx*ny/nThreads VECTOR CONTAINING THE HORIZONTAL 'i' INDICES OF THE SUBAPERTURES
	vector<vector<int>> vInd_j;			// [MT] nx*ny/nThreads VECTOR CONTAINING THE VERTICAL 'j' INDICES OF THE SUBAPERTURES
	vector<int>			vInd_i2;		// [MT] nx*ny - int(nx*ny/nThreads)*nThreads VECTOR CONTAINING THE HORIZONTAL 'i' INDICES OF THE REMAINING SUBAPERTURES
	vector<int>			vInd_j2;		// [MT] nx*ny - int(nx*ny/nThreads)*nThreads VECTOR CONTAINING THE VERTICAL 'j' INDICES OF THE REMAINING SUBAPERTURES
	vector<float>		vThresh;
	vector<int>			vEnabled;

	// CENTROID POSITIONS
	int				nDataSets;				// NUMBER OF nx*ny DATA CONTAINED IN THE SHARED MEMORY BUFFER 'SHM_Slopes'
	mcUINT32		local_counter;
	mcUINT32		display_image_count;
	vector<float>	vSH_cent_x;				// COMPUTED 'X' CENTROIDS FOR EACH SUBAPERTURE
	vector<float>	vSH_cent_y;				// COMPUTED 'Y' CENTROIDS FOR EACH SUBAPERTURE
	vector<float>	vSH_slope_x;			// COMPUTED 'X' SLOPES FOR EACH SUBAPERTURE
	vector<float>	vSH_slope_y;			// COMPUTED 'Y' SLOPES FOR EACH SUBAPERTURE
	vector<float>	vSH_phot;				// COMPUTED FLUX FOR EACH SUBAPERTURE
	vector<float>	vSH_ref_x;				// REFERENCE 'X' CENTROIDS FOR EACH SUBAPERTURE
	vector<float>	vSH_ref_y;				// REFERENCE 'Y' CENTROIDS FOR EACH SUBAPERTURE

	// ARRAY CONTAINING CENTROID POSITIONS FOR DISPLAY ONLY
	int				MaxModes;			// MAXIMUM NUMBER OF DM MODES TO CONSIDER (MUST BE A "MAX" BECAUSE OF UNRESIZABLE ASSOCIATED SHARED MEMORY BUFFER)
	//vector<float>	vFlatSHData;		// nx*ny*5*N VECTOR CONTAINING THE CENTROIDS, SLOPES AND FLUX HISTORY OF THE 'N' PROCESSED DATA BETWEEN EACH DISPLAY UPDATE (x FOLLOWED BY y FOLLWED BY FLUX)
	vector<float>	vModesDataHist;		// MaxModes*N VECTOR CONTAINING THE MODE PROJECTION HISTORY OF THE 'N' PROCESSED DATA BETWEEN EACH DISPLAY UPDATE (x FOLLOWED BY y FOLLWED BY FLUX)
	
	// TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	vector<vector<float>>	vCircBuf_FlatSHData;
	vector<vector<double>>	vCircBuf_Timestamps;
	vector<vector<float>>	vCircBuf_DM_cmd;
	size_t					nCircBuffers;
	size_t					circBufCounter;

	// DM MASK
	vector<int>		vSH_PupilMask;		// nx*ny VECTOR DEFINING VALID (1) OR INVALID (0) SUBAPERTURES
	vector<int>		vSH_PupilWeights;	// nx*ny VECTOR DEFINING 
	vector<int>		vDMPuptoSH;			// VECTOR CONTAINING THE DM SUBAPERTURE PROJECTION INDICES TO THE ACTUAL SH (nx*ny TOTAL SIZE)

	// PROJECTION MODE
	enum
	{
		AOC_PROJ_MODE_ZER = 1,
		AOC_PROJ_MODE_ZON
	} enumProjMode;
	int				projMode;

	// CALIBRATION MODES
	bool			bCalMode;
	int				nCalFrameCount;
	int				curCalMode;

	bool			bEnableCloseLoop;

	// DM COMMAND MODES AND DM OUTPUT RESPONSE MODES (SIZE DEPENDS OF THE NUMBER OF USED MODES)
	bool			bCalDataAvailable_ZER;	// FLAG SPECIFYING IF THE ZERNIKE CALIBRATION DATA ARE AVAILABLE
	int				nModes_ZER;				// NUMBER OF REGISTERED ZERNIKE MODES
	vector<float>	vDM_Cmd_ZER;			// VECTOR OF nModes*nAct ELEMENTS CONTAINING THE DM COMMANDS FOR EACH ZERNIKE MODE
	vector<float>	vDM_Modes_ZER_Slope;	// VECTOR OF nModes*nx*ny*2 ELEMENTS CONTAINING THE CORRESPONDING DM RESPONSE X FOLLOWED BY Y SLOPES
	vector<float>	vDM_Modes_ZER_iNorm;	// VECTOR OF nModes ELEMENTS CONTAINING THE SQUARED SUM OF EACH RECORDED SLOPE MODE (FOR NORMALIZATION OF PROJECTION ON EACH MODE)
	vector<float>	vDM_ProjOnModes_ZER;	// VECTOR OF nModes ELEMENTS CONTAINING THE RESULT OF THE PROJECTION OF THE CURRENT SLOPE MAP ONTO THE RECORDED MODES
	float			ZER_gain;
	//vector<double>	vZER_CmdMat;			// VECTOR CONTAINING THE ZERNIKE COMMAND MATRIX

	bool			bCalDataAvailable_ZON;	// FLAG SPECIFYING IF THE ZONAL CALIBRATION DATA ARE AVAILABLE
	int				nModes_ZON;				// NUMBER OF REGISTERED ZONAL MODES
	vector<float>	vDM_Cmd_ZON;			// VECTOR OF nModes*nAct ELEMENTS CONTAINING THE DM COMMANDS FOR EACH ZONAL ACTUATION
	vector<float>	vDM_Modes_ZON_Slope;	// VECTOR OF nModes*nx*ny*2 ELEMENTS CONTAINING THE CORRESPONDING DM RESPONSE X FOLLOWED BY Y SLOPES
	vector<float>	vDM_Modes_ZON_iNorm;	// VECTOR OF nModes ELEMENTS CONTAINING THE INVERSE (FOR FASTEST OPERATION) OF THE SQUARED SUM OF EACH RECORDED MODE (FOR NORMALIZATION OF PROJECTION ON EACH MODE)
	vector<float>	vDM_ProjOnModes_ZON;	// VECTOR OF nModes ELEMENTS CONTAINING THE RESULT OF THE PROJECTION OF THE CURRENT SLOPE MAP ONTO THE RECORDED MODES
	//vector<double>	vZON_CmdMat;			// VECTOR CONTAINING THE ZERNIKE COMMAND MATRIX

	vector<HANDLE>			vThreadHandles;
	vector<ThreadParams_CalcSHData>
							vTP_CalcSHData;
	vector<DWORD>			vThread_CalcSHDataID;

	// TIMESTAMPS
	vector<double>			timestamp;

	//LOG DIR INFO
	bool					isCal;

} ThreadParams_calcSH;

// FUNCTION PROTOTYPES
mcINT16 StartThread_calcSH(int _nThreads, int _DX, int _DY, int _nx, int _ny);

mcINT16 TerminateThread_calcSH();

//mcINT16 set_thread_params(int threadID, int _nsub, int *_v_xx, int *_v_yy, int *_vInd_i, int *_vInd_j, float *_SH_xtmp, float *_SH_ytmp, float *_SH_phot);
mcINT16 set_thread_params(int threadID, int _nsub, vector<int>& _v_xx, vector<int>& _v_yy, vector<int>& _vInd_i, vector<int>& _vInd_j, 
	vector<float>& _SH_xtmp, vector<float>& _SH_ytmp, vector<float>& _SH_slope_x, vector<float>& _SH_slope_y, vector<float>& _SH_phot, 
	vector<float>& _SH_refx, vector<float>& _SH_refy, vector<float>& _vThresh, vector<int>& _vEnabled);

int trigger_calc_SH_data_mt(long *vImg, float _bkg);
int trigger_proj_TT();
int trigger_proj_ZER();
int trigger_set_DM_correction_map();
int trigger_proj_ZON();
