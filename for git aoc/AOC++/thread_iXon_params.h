#pragma once

#include <cstdlib>
#include <Windows.h>

#include <MCODE.h>
#include <MCODE_SciLib.h>

using namespace MultiCODE;
using namespace std;

typedef struct
{
	// THREAD CREATION/TERMINATION EVENTS
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;

	// TRIGGERING EVENTS
	HANDLE			hShutdownCamera;
	HANDLE			hStartAcquisition;
	HANDLE			hStopAcquisition;
	HANDLE			hGetLastImage;
	
	// SIGNALING EVENTS
	HANDLE			hFrameReady;
	HANDLE			hAcquisitionDone;
	HANDLE			hExt_StatsUpdate;
	HANDLE			hExt_MissedFrame;
	HANDLE			hLiveStream;

	// SIMULATION MODE
	bool			bSimulation;

	// CAMERA CONFIGURATION PARAMETERS
	unsigned int	frame_count;
	double			time_diff_s;
	bool			bCamInitialized;
	bool			bStreamON;
	bool			bAcqON;
	bool			bCalibration;
	bool			bCoolingON;
	int				camID;
	int				startX;
	int				startY;
	int				DX;
	int				DY;
	int				binX;
	int				binY;
	float			expTime_s;
	float			accTime_s;
	float			kinTime_s;
	float			readoutTime;
	int				acqTimeout_ms;
	int				tempMin;
	int				tempMax;
	int				testLocalCount;

	// ACQUISITION LOOP WORING MODE (CALIBRATION, CLOSED LOOP, LIVE CAMERA)
	//bool			bCalMode;
	//bool			bClosedLoopMode;
	//int			nCalFrameCount;

	// LOW-LEVEL HARDWARE CONFIGURATION
	int				iVHSIndex;
	float			VSspeed;
	float			HZSspeed;
	int				readMode;
	int				acqMode;
	int				FTmode;						// FRAME TRANSFER MODE (0=OFF, 1=ON)
	int				trigMode;
	int				ampType;					// TYPE OF OUTPUT AMPLIFICATION (1: CONVENTIONAL, 0: EMCCD)
	int				ADCchannel;
	int				PreAmpGainIndex;

	// IMAGE PROCESSING VARIABLES
	vector<long>	vImgBufferCopy;
	bool			bAsyncProc;

	// THREAD VARIABLES
	bool			bEndThread;
	DWORD			WaitTimeout;

	vector<long>				vImgBuffer;

	//DWORD			curFrame;
	DWORD			lastFrame;
	DWORD			internalFrameCount;
	LARGE_INTEGER	li_timestamp;
	double			FPS;
	//int				FPSinterval;					// NUMBER OF FRAME FOR COMPUTING THE AVERAGE
	double			display_image_interval_s;		// TIME INTERVAL BETWEEN TWO DISPLAY BUFFER UPDATES
	int				display_image_count;			// NUMBER OF IMAGES BETWEEN TWO DATA MONITORINGS (COMPUTED FROM display_image_interval_s)
	int				nCalFrames;

	MCODE_SciLib	sci;

	// CRITICAL SECTIONS
	CRITICAL_SECTION	CSCopyImgBuffer;
} ThreadParams_iXon;

typedef struct {
	// THREAD CREATION/TERMINATION EVENTS
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;

	// TRIGGERING EVENTS
	HANDLE			hStartLiveStream;
	HANDLE			hStopLiveStream;

	// SIGNALING EVENTS
	HANDLE			hLiveStreamDone;
	HANDLE			hMissedFrame;

	// THREAD VARIABLES
	bool bEndThread;
	DWORD WaitTimeout;
	int maxImages;
	int currentIndex;
	unsigned int imgCount;
	bool isMissed;

	string fitsDir;
} ThreadParams_Simu_iXon;

// CLASS FOR LOADING FITS FILES

class SimulationLoader {
public:
	static SimulationLoader& getInstance() {
		static SimulationLoader simulationLoader;
		return simulationLoader;
	}

	void loadFullRepo(const string &repoPath) {
		DIR* dir;
		struct dirent* ent;
		if ((dir = opendir(repoPath.c_str())) != NULL) {
			while ((ent = readdir(dir)) != NULL) {
				string filename(ent->d_name);
				if (filename.substr(filename.find_last_of(".") + 1) == "fits") {
					add(repoPath + filename);
				}
			}
			closedir(dir);
		}
		else {
			cout << "could not open repo sor simulation loader" << endl;
		}
	}

	void add(const string& filepath) {
		MCODE_Export_FITS FITS;

		vector<long> vBuffer;
		vector<int> vSize;

		string test = filepath;
		if (FITS.Load(test, vBuffer, vSize) != SUCCESS) {
			cout << "Error when loading fits in simulation loader" << endl;
		}
		else {
			cout << test << " successfully load" << endl;
		}
		_values.push_back(vBuffer);
	}

	vector<long>& get(int index) {
		return _values.at(index);
	}

	vector<vector<long>>& getAllValues() {
		return _values;
	}

	void emptyValues() {
		_values.clear();
	}

	void operator=(const SimulationLoader&) = delete;
private:
	SimulationLoader() {};
	vector<vector<long>> _values;
};

// END CLASS FOR LOADING FITS FILE

mcINT16 AOC_SetHWConf(ThreadParams_iXon *param);
mcINT16 AOC_InitializeCamera(ThreadParams_iXon *param);
mcINT16 AOC_ShutdownCamera(ThreadParams_iXon *param);
mcINT16 AOC_ProcessSHData(ThreadParams_iXon *param);

mcINT16 GetCS_CopyImgBuffer(CRITICAL_SECTION *_CS);

mcINT16 StartThread_iXon();
mcINT16 TerminateThread_iXon();
static mcUINT32 Thread_iXon(void *pParam);

mcINT16 StartThread_Simu_iXon();
mcINT16 TerminateThread_Simu_iXon();
static mcUINT32 Thread_Simu_iXon(void* pParam);