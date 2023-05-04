/* =========================================================================
   Source code for the camera server running the andor iXon camera for CIAO.
   The server awaits from simple commands like "stream", "quit", "tint xxx"
   on a named pipe (aka a fifo, for "first in, first out") and executes them
   assuming that the command is valid.
   
   Implement the revised shared memory data structure developed during
   the summer 2017.


   Frantz Martinache.
   ========================================================================= */

#include <stdio.h>
#include <AclAPI.h>

#ifdef __GNUC__
#  if(__GNUC__ > 3 || __GNUC__ ==3)
#	define _GNUC3_
#  endif
#endif

#define WINDOWS
#define _GNUC3_

#if defined(_GNUC3_)
#  include <iostream>
#  include <fstream>
#include <io.h>
   using namespace std;
#else
#  include <iostream.h>
#  include <fstream.h>
#endif

#include <string.h>
#include <stdlib.h>
#if defined(LINUX)
	#include <unistd.h>
	#include <sys/mman.h>
	#include "atmcdLXd.h"
#elif defined(WINDOWS)
	#include <atmcd32d.h>
#endif
#include <sys/stat.h>
#include <fcntl.h>

#include <fitsio.h>

//#define _TIMESPEC_DEFINED
//#include "SHM_DS.h"
#undef _TIMESPEC_DEFINED
#include <pthread.h>

#include "ImageStruct.h"
#include "ImageCreate.h"

#include <time.h>
#include <process.h>

#pragma comment(lib,"Advapi32.lib")

// ==========================================================================
//                              CONSTANTS 
// ==========================================================================
#define LINESZ 256                         // string length
#if defined(LINUX)
	#define IXON_SM_FNAME "/tmp/ixon.shm"      // shared memory file for iXon
#else
	#define IXON_SM_FNAME "R:\\ixon.shm"      // shared memory file for iXon
#endif
#define msleep(x) usleep((int)(1000*(x))); // sleep in milli-seconds

// ==========================================================================
//                  structure storing acquisition information
// ==========================================================================
typedef struct {
  float exp_time; // exposure time in seconds
  float acc_time; // accumulate cycle time in seconds (cf SDK)
  float kin_time; // kinetic cycle time in seconds (cf SDK)
  // ---------
  int XW;         // image (window) width
  int YW;         // image (window) height
  // ---------
  int nleft;      // number of images left in acquisition
  int nextframe;  // index of the next frame to write
  // ---------
  bool bstreamon; // camera streaming? (continuous acq.)
  bool bacqon;    // camera acquiring? 
  bool babort;    // abort command was issued?
  bool bcamOK;    // the happy camera flag!
} cam_config;



// THREAD PARAMETERS DEFINITIONS
#define THREAD_CREATION_TIMEOUT		5000

typedef struct _ThreadParams_Acquisition
{
	// CREATION/TERMINATION EVENTS
	HANDLE		hThreadCreated;
	HANDLE		hTerminateThread;
	
	// TRIGGERING EVENTS
	HANDLE		hStopAcquisition;
	HANDLE		hStartAcquisition;

	// SIGNALING EVENTS
	HANDLE		hAcquisitionDone;
	HANDLE		hAcquisitionStopped;

	float		ExpTime_s;

	bool		bEndThread;

	cam_config	*pCamConfig;
} ThreadParams_Acquisition;

HANDLE hThreadAcquisition;
DWORD ThreadAcquisitionID;

HANDLE hAndor_statusEvent;

int	StartThread_Acquisition(ThreadParams_Acquisition *param);
UINT32 Thread_Acquisition(void *param);

// ==========================================================================
//                          GLOBAL VARIABLES
// ==========================================================================
cam_config *cconf;                       // pointer to acq. info structure

char fits_name[LINESZ] = "test.fits";    // temporary feature for tests
char conf_name[LINESZ] = "config.xml";   // OptAcqr config file

#if defined(LINUX)
char myfifoin[LINESZ] = "/home/ciaodev/bin/ixon_fifo_in";  // out-going pipe
char myfifout[LINESZ] = "/home/ciaodev/bin/ixon_fifo_out"; // out-going pipe
#else
// THE PIPE NAME MUST BE IN THE FORM \\.\pipe\pipename (https://docs.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-createnamedpipea)
char myfifoin[LINESZ] = "\\\\.\\pipe\\ixon_fifo_in";  // out-going pipe
char myfifout[LINESZ] = "\\\\.\\pipe\\ixon_fifo_out"; // out-going pipe
#endif

//int first_im = 0;
//int last_im = 0;

// size_t shm_size = 0;                     // shared memory size in bytes
// int fd_shm = -1;                         // file descriptor for shared memory
// SHM_DATA *shm_data;                      // pointer to shared memory for iXon

IMAGE *imarray;            // pointer to image

// ==========================================================================
//                              PROTOTYPES
// ==========================================================================
int CameraSelect(int iNumArgs, char* szArgList[]);
void cam_config_init(cam_config *camconf);
void printBits(unsigned int num);

// ==========================================================================
//                         FUNCTIONS 
// ==========================================================================
void printBits(unsigned int num) {
  for(unsigned int bit=0;bit<(sizeof(unsigned int) * 8); bit++) {
    printf("%i ", num & 0x01);
    num = num >> 1;
  }
  printf("\n");
}

// ---------------------------------------------------------
//        initialize the camconf structure
// ---------------------------------------------------------
void cam_config_init(cam_config* camconf) {
  camconf->exp_time = 0.00001; // default: shortest exp time
  camconf->kin_time = 0.0;
  camconf->bstreamon = false;
  camconf->bacqon = false;
  camconf->babort = false;
  camconf->bcamOK = true;
}

// ==========================================================================
//                         ACQUISITION THREAD
// ==========================================================================
void* acquire(void *params)
{
	// camera acquisition thread
	cam_config* camconf = (cam_config*) params;
	int status = 0;
	int nel = camconf->XW * camconf->YW;
	int nfr = 100, ifr = 0;     // variables used to estimate frame rate
	float t0 = 0.0, t1 = 0.0;  // time variables for frame rate
	int idisp = 0, ndisp = 100; // control the text output refresh rate
	double ifrate = 0.;              // frame rate (average over nfr frames)
	double itiming = 0.;
	struct mytimespec now;       // clock readout
	struct tm* ptm;
	float *timing = (float*) malloc(nfr * sizeof(float));

#if defined(LINUX)
	int first_im, last_im;       // indices for avail. images in circ buffer
	int last_read_im;            // index of last read image
	int first_valid, last_valid; // indices of valid images
#elif defined(WINDOWS)
	long first_im;				// indices for avail. images in circ buffer
	long last_im;
	long cnt_first;
	long last_read_im;            // index of last read image
	long first_valid, last_valid; // indices of valid images
#endif
	unsigned long error;         // error code

	last_read_im = -1;           // keep track of last read image index
	for (ifr = 0; ifr < nfr; ifr++) timing[ifr] = 0.1; // init timing array
	ifr = 1;

	// acquisition loop
	printf("\n");

	if (camconf->bstreamon)
	{
		StartAcquisition();
	}
 
//      // what is the time ??
//#if defined(WINDOWS)
//    myclock_gettime(0, &now);
//	//const time_t now_sec = now.tv_sec;
//	time_t now_sec = now.tv_sec;
//	ptm = std::gmtime(&now_sec);
//#elif defined(LINUX)
//    clock_gettime(CLOCK_REALTIME, &now);
//	ptm = gmtime(&(now.tv_sec));
//#endif
//
//  t0 = (float)(now.tv_nsec)*1e-9 + (float)(ptm->tm_sec);

	LARGE_INTEGER li_freq;
	LARGE_INTEGER li_t0;
	LARGE_INTEGER li_now;
	double pc_freq = 0.;
	double pc_start = 0.;

	QueryPerformanceFrequency(&li_freq);
	double iFreq = 1./double(li_freq.QuadPart);
	double avg_delta = 0.;
	double meas_frate = 0.;
	int counter = 0;

	long delta = 0;
	DWORD camStatus = 0;

	QueryPerformanceCounter(&li_t0);

	// TAKE A FIRST IMAGE TO RECORD THE COUNTER
	GetNumberNewImages(&first_im, &last_im);
	cnt_first = last_im;

	while (camconf->nleft > 0)
	{
		// CODE FOR 'Internal trigger mode' MODE
		error = GetNumberNewImages(&first_im, &last_im);
		//WaitForSingleObject(,100);
		//GetStatus(&status);
		//printf("Status: %d\n",status);
		//WaitForAcquisition();

		//if (true)
		//{
		//	switch(WaitForSingleObject(hAndor_statusEvent,1000))
		//	{
		//	case WAIT_OBJECT_0:
		//		GetCameraEventStatus (&camStatus);
		//		break;
		//	case WAIT_TIMEOUT:
		//		printf("TIMEOUT!!\n");
		//		break;
		//	}
		//}

		//if (camStatus == 1)
		//{
		//	QueryPerformanceCounter(&li_now);
		//	avg_delta += double(li_now.QuadPart - li_t0.QuadPart) * iFreq;
		//	li_t0 = li_now;
		//	counter++;
		//}

		//if (counter == 100)
		//{
		//	printf("Timing: %d %.3e\n", counter, avg_delta/double(counter));
		//	counter = 0;
		//	avg_delta = 0.;
		//}

		//// CODE FOR 'Software trigger' MODE
		//SendSoftwareTrigger();
		//WaitForAcquisition();
		//GetNumberNewImages(&first_im, &last_im);

		if (last_im > last_read_im)
		{
			// a new image is available!
			delta += last_im - last_read_im;
			//printf("delta image: %d %d %d %d\n",last_im - last_read_im, last_im, last_read_im,error);
			last_read_im = last_im;

			//      // what is the time ??
			//#if defined(WINDOWS)
			//    myclock_gettime(0, &now);
			//	const time_t now_sec = now.tv_sec;
			//	ptm = std::gmtime(&now_sec);
			//#elif defined(LINUX)
			//    clock_gettime(CLOCK_REALTIME, &now);
			//	ptm = gmtime(&(now.tv_sec));
			//#endif
			//
			//      t1 = (float)(now.tv_nsec)*1e-9 + (float)(ptm->tm_sec);

//			QueryPerformanceCounter(&li_now);

			// estimate the frame rate
//			timing[ifr] = double(li_now.QuadPart - li_t0.QuadPart) * iFreq;
//			li_t0 = li_now;
			//ifr++;
			//if (ifr == nfr) ifr = 0;
      
			//ifrate = 0.0;
			//for (int i = 0; i < nfr; i++)
			//{
			//	ifrate += timing[i];
			//}

			//ifrate = (float)(nfr) / ifrate;


			// write to shared memory
			imarray[0].md[0].write = 1;                  // signal you are about to write
			//error = GetImages(last_im, last_im, (long*)imarray[0].array.SI32, nel, &first_valid, &last_valid);
			//error = GetNewData((long*)imarray[0].array.SI32, nel);
			error = GetMostRecentImage((long*)imarray[0].array.SI32, nel);
			//Sleep(1);
			imarray[0].md[0].write = 0;                  // signal done writing data
			imarray[0].md[0].cnt0 ++;                    // increment counter!
			imarray[0].md[0].cnt1 ++;                    // increment counter!

			idisp++;
		}

		// check for abort
		if (camconf->babort)
		{
			camconf->nleft = 0;
			camconf->bstreamon = false;
			camconf->babort = false;
		}

		if (idisp == ndisp)
		{
			QueryPerformanceCounter(&li_now);
			ifrate = (last_im - cnt_first) / (double(li_now.QuadPart - li_t0.QuadPart) * iFreq);
			itiming = 1./ifrate;
			avg_delta = double(delta) / double(idisp);
			meas_frate = 1. / (avg_delta * itiming);
			printf("%02d Image [%9ld] (i0=%6d,i1=%6d) int.tmg=%.3e iFPS=%.1fHz avg.delta=%.3e mFPS=%.1f, err=%ld\n", 
				idisp, imarray[0].md[0].cnt0, first_im, last_im, itiming, ifrate, avg_delta, meas_frate, error);
			//fflush(stdout);
			cnt_first = last_im;
			idisp = 0;
			li_t0 = li_now;
			delta = 0;
		}

		//msleep(1); // sleep in milli-seconds
		if (!camconf->bstreamon) camconf->nleft--; // decrement if NOT streaming

	}

	//t1 = (float)(now.tv_nsec)*1e-9 + (float)(ptm->tm_sec);
	//printf("\n");
	//double avg_time = (t1 - t0)/ifr;
	//double avg_freq = 1./avg_time;
	//printf("Average time between images (freq): %.3e (%.3f)\n",avg_time,avg_freq);

	camconf->bacqon    = false; // updating acquisition flags
	camconf->bstreamon = false; // updating acquisition flags
	return NULL;
}

// ==========================================================================
//                              MAIN PROGRAM
// ==========================================================================
int main(int argc, char* argv[]) {
	pthread_t tid_acqr;      // thread id for acquisition

	int mytest = -1;         // useful for ... tests
	long naxis = 2;          // number of axes
	uint8_t atype;           // data type code
	uint32_t *imsize;        // image size
	int shared;              // set to 1 if image is in shared memory
	int NBkw;                // number of keywords supported

	AndorCapabilities *caps; // TEST

	unsigned long error;
	bool quit;
	float fChoice;
	int iChoice;
	int gain = -1;
	int temp;
	int temp_min = 1, temp_max = 1;
	unsigned int state;

	int glow = -1, ghigh = -1;
	int nb_pre_gains = -1;
	int nb_adcs = -1, nb_amps = -1, nb_speeds = -1, 
	adc_channel = -1, amp_type = -1;
	float pre_gain = 0.0, hzs_speed = 0.0;

	char *cmd_in = (char*) malloc(LINESZ * sizeof(char));
	char *buf = (char*) malloc(LINESZ * sizeof(char));
	char *auxstr = (char*) malloc(LINESZ * sizeof(char));
	char *msgout = (char*) malloc(LINESZ * sizeof(char));
	mytest *= 1;

	cconf = (cam_config*) malloc(sizeof(cam_config));
	cam_config_init(cconf);

	//sem_t *test_sem = (sem_t*)my_sem_open("test_semaphore",O_CREAT,0644,1);


  // ------------------------------------------
  //            Initialize CCD
  // ------------------------------------------
	if (CameraSelect (argc, argv) < 0) {
		cout << "*** CAMERA SELECTION ERROR" << endl;
		return -1;
	}
  
#if defined(WINDOWS)
	error = Initialize("");
#elif defined(LINUX)
	error = Initialize((char*)"/usr/local/etc/andor/");
#endif
	if(error != DRV_SUCCESS){
		cout << "Initialisation error...exiting" << endl;
		return(1);
	}
  
  // ------------------------------------------
  //    startup configuration of the camera
  // ------------------------------------------
#if defined(LINUX)
	sleep(2);                         // sleep to allow initialization to complete
#elif defined(WINDOWS)
	Sleep(2000);
#endif
	error = SetPCIMode(1,1);
	hAndor_statusEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
	////error = SetDriverEvent(hAndor_statusEvent);
	error = SetAcqStatusEvent(hAndor_statusEvent);
	error = SetDMAParameters(1, 0.001);
	SetFrameTransferMode(1);			// Frame transfer mode is OFF
	SetReadMode(4);                   // Set Read Mode to --Image--
	SetAcquisitionMode(5);            // Set Acq. mode to --Run till abort--
	SetShutter(1,1,50,50);            // Initialize Shutter
	SetExposureTime(cconf->exp_time);           // Set initial exposure time
	GetDetector(&(cconf->XW), &(cconf->YW));    // Get Detector dimensions  
	SetImage(1, 1, 1, cconf->XW, 1, cconf->YW); // Setup Image dimensions

	SetKineticCycleTime(cconf->kin_time); // test
	error = SetNumberAccumulations(1);

	error = GetAcquisitionTimings(&(cconf->exp_time), 
				&(cconf->acc_time), 
				&(cconf->kin_time));
	printf("\n >> Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s\n",
		cconf->exp_time, cconf->acc_time, cconf->kin_time);

	caps = (AndorCapabilities *) malloc(sizeof(AndorCapabilities));
	caps->ulSize = sizeof(AndorCapabilities);
	error = GetCapabilities(caps);
	if(error != DRV_SUCCESS)
		cout << "GET Capabilities fail" << endl;

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
		error = GetPreAmpGain(i, &pre_gain);
		printf("%d: gain = %.2f\n", i, pre_gain);
	}
	error = SetPreAmpGain(0); // set pre-amp gain to 1.0

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

	amp_type    = 0; // type of output amplification (1: conventional, 0: EMCCD)
	adc_channel = 0; //

	error = SetOutputAmplifier(amp_type);
	error = SetADChannel(adc_channel);

	printf("For AD channel = %d: ", adc_channel);
	error = GetNumberHSSpeeds(adc_channel, amp_type, &nb_speeds);
	printf("%d speeds are available in this mode\n", nb_speeds);

	for (int i = 0; i < nb_speeds; i++)
	{
		error = GetHSSpeed(adc_channel, amp_type, i, &hzs_speed);
		printf("- i=%d: speed = %.1f MHz\n", i, hzs_speed);
	}
	printf("\n");

	// ====================================================================

	//SetHSSpeed(amp_type, 0);
	//  SetOutputAmplifier(amp_type);

	//SetHighCapacity(1); // 
	error = GetEMCCDGain(&gain);

	if(error != DRV_SUCCESS)
	cout << "GET EMCCDGain fail" << endl;

	printf("EMCCD Gain = %d\n", gain);

	GetEMGainRange(&glow, &ghigh);
	printf("EM Gain range = %d - %d\n", glow, ghigh);

	GetMCPGain(&gain);
	printf("MCP   Gain = %d\n", gain);

	GetMCPGainRange(&glow, &ghigh);
	printf("MCP Gain range = %d - %d\n", glow, ghigh);

	cconf->exp_time = 0.00001;
	cconf->kin_time = 0.0;
	SetExposureTime(cconf->exp_time);
	SetKineticCycleTime(cconf->kin_time);
	float ReadoutTime = -1.f;
	GetReadOutTime(&ReadoutTime);

	int trigMode = 0;				// 0: Internal trigger | 10: Software trigger
	SetTriggerMode(trigMode);
	if (IsTriggerModeAvailable(trigMode) != DRV_SUCCESS)
	{
		printf("Trigger mode %d is not available!\n",trigMode);
	}

	if (caps->ulTriggerModes & AC_TRIGGERMODE_CONTINUOUS)
	{
		printf("Correct Hardware for Continuous Mode.\n");
	}
	else
	{
		printf("InCorrect Hardware for Continuous Mode.\n");
	}

	float fs;
	int iVHSIndex = 1;
	//GetFastestRecommendedVSSpeed (&iVHSIndex,&fs);
	SetVSSpeed(iVHSIndex);

	float speed = -1.f;
	GetVSSpeed(iVHSIndex, &speed);
	printf("Vertical Speed set to %g microseconds per pixel shift\r\n",speed);

	if (caps->ulSetFunctions & AC_SETFUNCTION_BASELINECLAMP)
	{
		error = SetBaselineClamp(1);
		if(error != DRV_SUCCESS)
		{
			printf("Set Baseline Clamp Error\n");
		}
	}


	// ------------------------------------------
	//    special setup for the temperature
	// ------------------------------------------
	if (false)
	{
		state = CoolerON();
		cout << "CoolerON status code: " << state << endl;

		state = GetTemperatureRange(&temp_min, &temp_max);
		cout << "Detector temperature range = " << temp_min << " - " << temp_max << endl;

		state = GetTemperature(&temp);    // get detector temperature
		cout << "Detector temperature = " << temp << " deg" << endl;

		temp = 0; // temperature setting (air cooling)
		state = SetTemperature(temp);

		state = GetTemperature(&temp);    // get detector temperature
		cout << "Detector temperature = " << temp << " deg" << endl;

		cout << "GetTemperature status code: " << state << endl;

		while(state != DRV_TEMPERATURE_STABILIZED)
		{
			Sleep(1000);
			state = GetTemperature(&temp);    // get detector temperature
			printf("\rTemperature = %+03d deg - status code: %6d", temp, state);
			fflush(stdout);
		}

		switch(state)
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
		state = CoolerOFF();
		cout << "CoolerOFF status code: " << state << endl;
		cout << "Current temperature: " << fTemp << endl;
	}

	// ------------------------------------------------
	// create pipes for interaction with other programs
	// ------------------------------------------------
#if defined(LINUX)
	if (mkfifo(myfifoin, 0777) != 0) printf("Could not create fifo!\n");
	if (mkfifo(myfifout, 0777) != 0) printf("Could not create fifo!\n");
#elif defined(WINDOWS)
	HANDLE hPipeIn;
	HANDLE hPipeOut;

	DWORD PipeWait = 1; // NMPWAIT_USE_DEFAULT_WAIT OR NUMBER OF MILLISECONDS FOR TIMEOUT

	DWORD dwWritten;
	hPipeIn = CreateNamedPipe(TEXT(myfifoin),
		PIPE_ACCESS_DUPLEX,
		PIPE_TYPE_MESSAGE | PIPE_READMODE_BYTE | PIPE_WAIT,   // FILE_FLAG_FIRST_PIPE_INSTANCE is not needed but forces CreateNamedPipe(..) to fail if the pipe already exists...
		1,
		1024 * 64,
		1024 * 64,
		PipeWait,
		NULL);

	hPipeOut = CreateNamedPipe(TEXT(myfifout),
		PIPE_ACCESS_DUPLEX,
		PIPE_TYPE_MESSAGE | PIPE_READMODE_BYTE | PIPE_WAIT,   // FILE_FLAG_FIRST_PIPE_INSTANCE is not needed but forces CreateNamedPipe(..) to fail if the pipe already exists...
		1,
		1024 * 64,
		1024 * 64,
		PipeWait,
		NULL);


#endif

	// ------------------------------------------
	//          setup shared memory 
	// ------------------------------------------
	naxis     = 2;                // 2D image
	imarray   = (IMAGE*) malloc(sizeof(IMAGE));
	imsize    = (uint32_t *) malloc(naxis * sizeof(uint32_t));
	imsize[0] = cconf->XW;       // image dimensions
	imsize[1] = cconf->YW;       // image dimensions
	atype     = _DATATYPE_INT32; // camera SDK writes "int"
	shared    = 1;               // image will be in shared memory
	NBkw      = 10;              // allocate space for 10 keywords
	ImageCreate(&imarray[0], "ixon", naxis, imsize, atype, shared, NBkw);


	// ------------------------------------------
	//          main server loop
	// ------------------------------------------
	int rfd = 0; // file descriptor for input fifo
	int wfd = 0; // file descriptor for output fifo
	int foo = -1;

	ThreadParams_Acquisition TP_ThreadAcquisition;
	TP_ThreadAcquisition.pCamConfig = cconf;
	
	if (StartThread_Acquisition(&TP_ThreadAcquisition) == 0)
	{
		printf("Thread_Acquisition started SUCCESSFULLY.\n");
		quit = false;
	}
	else
	{
		printf("Thread_Acquisition FAILED to start.\n");
		quit = true;
	}

	quit = false;

	printf("iXon server ready.\n");
	fflush(stdout);

#if defined(LINUX)
	do
	{
	// ==============================================
	//    listen for commands on the input fifo
	// ==============================================
	rfd = open(myfifoin, O_RDONLY);
	sprintf(cmd_in, "%s", "          ");
	read(rfd, cmd_in, LINESZ);
	close(rfd);
	rfd = 0;

	// ===============================================
	//       processing the command from pipe
	// ===============================================

	// -----------------------------------
	//      abort current acquisition
	// -----------------------------------
	if (strncmp(cmd_in, "abort", 5) == 0)
	{
		if (cconf->bacqon)
		{
			cconf->babort = true;
			AbortAcquisition(); // camera command
		}
	}

	// -----------------------------------
	//         close the program
	// -----------------------------------
	if (strncmp(cmd_in, "quit", 4) == 0) {
		if (cconf->bacqon)
		{
			cconf->babort = true;
			AbortAcquisition(); // camera command
			printf("\nQuitting program\n");
#if defined(LINUX)
			sleep(5);
#elif defined(WINDOWS)
			Sleep(5000);
#endif
		}
		quit = true;
	}

	// -----------------------------------
	//       take a single image
	// -----------------------------------
	if (strncmp(cmd_in, "take ", 5) == 0) {
		foo = -1;
		sscanf(cmd_in, "%s %d", auxstr, &foo);
		cconf->nleft = (foo < 0) ? 1: foo; 
		cconf->bacqon = true;
		pthread_create(&tid_acqr, NULL, acquire, cconf);
	}

	// -----------------------------------
	//       continuous acquisition
	// -----------------------------------
	if (strncmp(cmd_in, "stream", 6) == 0)
	{
		if (!cconf->bstreamon)
		{
			cconf->bacqon = true;
			cconf->bstreamon = true;
			cconf->nleft = 1;
			pthread_create(&tid_acqr, NULL, acquire, cconf);
		}
	}

    // -----------------------------------
    //       update exposure time
    // -----------------------------------
    if (strncmp(cmd_in, "tint ", 5) == 0) {
		sscanf(cmd_in, "%s %f", auxstr, &fChoice);
		//printf("Exposure time: %f\n", fChoice);
		AbortAcquisition(); // camera command

		cconf->exp_time = fChoice;
		cconf->kin_time = 0.0;
		SetExposureTime(cconf->exp_time);
		SetKineticCycleTime(cconf->kin_time);

		error = GetAcquisitionTimings(&(cconf->exp_time), 
			&(cconf->acc_time), 
			&(cconf->kin_time));

		printf("Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s\n",
			cconf->exp_time, cconf->acc_time, cconf->kin_time);

#if defined(LINUX)
		sleep(1);
#elif defined(WINDOWS)
		Sleep(1000);
#endif
    }

    // -----------------------------------
    //       query exposure time
    // -----------------------------------
	if (strncmp(cmd_in, "emgain ", 7) == 0) {
		sscanf(cmd_in, "%s %d", auxstr, &iChoice);
		AbortAcquisition(); // camera command

		error = SetEMCCDGain(iChoice);
#if defined(LINUX)
		sleep(1);
#elif defined(WINDOWS)
		Sleep(1000);
#endif
	}

	// -----------------------------------
	//       query exposure time
	// -----------------------------------
	if (strncmp(cmd_in, "tint?", 5) == 0)
	{
		wfd = open(myfifout, O_WRONLY);
		sprintf(msgout, "%.5f", cconf->exp_time);
		write(wfd, msgout, LINESZ);
		close(wfd);
		wfd = 0;
	}

  } while(!quit);

#elif defined(WINDOWS)
	DWORD dwRead;

	//hPipeIn = INVALID_HANDLE_VALUE;

	if (hPipeIn != INVALID_HANDLE_VALUE)
	{

		while((!quit) || (ConnectNamedPipe(hPipeIn, NULL) != FALSE))   // wait for someone to connect to the pipe
		{
			//DWORD dw = GetLastError(); 
			//printf("failed with error %d\n",dw);
			while (!quit)
			{
				if (ReadFile(hPipeIn, cmd_in, sizeof(cmd_in) - 1, &dwRead, NULL) != FALSE)
				{
					cmd_in[dwRead] = '\0';

					int len = dwRead;
					
					DWORD bytesToRead;
					DWORD bytesAvail;
					DWORD bytesLeft;
					// THIS METHOD CHECKS WHETHER A PIPE STILL HAS BYTES TO READ FROM THE ONGOING MESSAGE
					PeekNamedPipe(hPipeIn, &buf, sizeof(buf),&bytesToRead,&bytesAvail,&bytesLeft);
					while(bytesToRead != 0)
					{
						ReadFile(hPipeIn, cmd_in + len, sizeof(cmd_in) - 1, &dwRead, NULL);
						len += dwRead;
						PeekNamedPipe(hPipeIn, &buf, sizeof(buf),&bytesToRead,&bytesAvail,&bytesLeft);
					}

					/* add terminating zero */
					cmd_in[len] = '\0';

					/* do something with data in buffer */
					// ===============================================
					//       processing the command from pipe
					// ===============================================

					// -----------------------------------
					//      abort current acquisition
					// -----------------------------------
					if (strncmp(cmd_in, "abort", 5) == 0)
					{
						if (cconf->bacqon)
						{
							cconf->babort = true;
							if (WaitForSingleObject(TP_ThreadAcquisition.hAcquisitionDone,5000) == WAIT_TIMEOUT)
							{
								printf("TIMEOUT while waiting for acquisition routine to terminate\n");
							}

							AbortAcquisition(); // camera command
						}
					}

					// -----------------------------------
					//         close the program
					// -----------------------------------
					if (strncmp(cmd_in, "quit", 4) == 0)
					{
						if (cconf->bacqon)
						{
							cconf->babort = true;
							if (WaitForSingleObject(TP_ThreadAcquisition.hAcquisitionDone,5000) == WAIT_TIMEOUT)
							{
								printf("TIMEOUT while waiting for acquisition routine to terminate\n");
							}
							AbortAcquisition(); // camera command
							printf("\nQuitting program\n");
							Sleep(5000);
						} 
						quit = true;
					}

					// -----------------------------------
					//       set temperature goal
					// -----------------------------------
					if (strncmp(cmd_in, "tempg ", 6) == 0)
					{
						sscanf(cmd_in, "%s %d", auxstr, &iChoice);
						if ((iChoice < -60) || (iChoice > 20))
						{
							printf("\nTemperature set point (%d) not in range (-60<T<+20). Nothing done.\n",iChoice);
						}
						else
						{
							SetTemperature(iChoice); // camera command
							printf("\nTempertaure set point set to %d Cdeg.\n",iChoice);
						}
						Sleep(1000);
					}


					// -----------------------------------
					//       set temperature goal
					// -----------------------------------
					if (strncmp(cmd_in, "cooling ", 8) == 0)
					{
						sscanf(cmd_in, "%s %d", auxstr, &iChoice);
						if (iChoice == 0)
						{
							state = CoolerOFF();
							printf("\nCooling is now OFF (error=%d).\n",state);
						}
						else if (iChoice == 1)
						{
							state = CoolerON();
							printf("\nCooling is now ON (error=%d).\n",state);
						}
						else
						{
							printf("\nUnhandled 'cooling' value (value=%d, valid=0|1).\n",iChoice);
						}
						Sleep(1000);
					}


					if (strncmp(cmd_in, "fTemp?", 6) == 0)
					{
						bool flag = true;
						float fTemp = -999.f;
						state = GetTemperatureF(&fTemp);    // get detector temperature

						switch(state)
						{
						case DRV_ACQUIRING:
							cout << "Not available (acquisition in progress)." << endl;
							fTemp = -999.f;
							flag = false;
							break;
						case DRV_TEMP_STABILIZED: 
							cout << "STABILIZED." << endl;
							break;
						case DRV_ERROR_ACK:
							cout << "ERROR: Unable to communicate with card." << endl;
							fTemp = -999.f;
							flag = false;
							break;
						case DRV_TEMP_DRIFT:
							cout << "Temperature DRIFT." << endl;
							break;
						case DRV_TEMP_NOT_STABILIZED:
							cout << "NOT STABILIZED." << endl;
							break;
						case DRV_TEMP_OFF: 
							cout << "Cooler OFF." << endl;
							break;
						case DRV_TEMP_NOT_REACHED: 
							cout << "COOLING in progress..." << endl;
							break;
						default:
							cout << "Cooler status unknown." << endl;
						}

						if (flag == true)
						{
							cout << "Detector temperature = " << fTemp << " deg" << endl;
						}

						std::sprintf(msgout, "%.3f", fTemp);
						if (WriteFile(hPipeIn,
							msgout,
							std::strlen(msgout),   // = length of string + terminating '\0' !!!
							&dwWritten,
							NULL) == FALSE)
						{
							printf("ERROR: WriteFile(hPipeIn) failed!\n");
						}
						FlushFileBuffers(hPipeIn);
	
						Sleep(1000);
					}


					// -----------------------------------
					//       take a single image
					// -----------------------------------
					if (strncmp(cmd_in, "take ", 5) == 0)
					{
						foo = -1;
						sscanf(cmd_in, "%s %d", auxstr, &foo);
						cconf->nleft = (foo < 0) ? 1: foo; 
						cconf->bacqon = true;
						pthread_create(&tid_acqr, NULL, acquire, cconf);
					}

					// -----------------------------------
					//       continuous acquisition
					// -----------------------------------
					if (strncmp(cmd_in, "stream", 6) == 0)
					{
						if (!cconf->bstreamon)
						{
							cconf->bacqon = true;
							cconf->bstreamon = true;
							cconf->nleft = 1;
							//pthread_create(&tid_acqr, NULL, acquire, cconf);
							SetEvent(TP_ThreadAcquisition.hStartAcquisition);
						}
					}

					// -----------------------------------
					//       update exposure time
					// -----------------------------------
					if (strncmp(cmd_in, "tint ", 5) == 0)
					{
						sscanf(cmd_in, "%s %f", auxstr, &fChoice);
						//printf("Exposure time: %f\n", fChoice);
						AbortAcquisition(); // camera command

						cconf->exp_time = fChoice;
						cconf->kin_time = 0.0;
						SetFrameTransferMode(0);				// Frame transfer mode is OFF
						SetReadMode(4);							// Set Read Mode to --Image--
						SetAcquisitionMode(5);					// Set Acq. mode to --Run till abort--
						SetExposureTime(cconf->exp_time);
						SetKineticCycleTime(0.f);

						error = GetAcquisitionTimings(&(cconf->exp_time), 
									&(cconf->acc_time), 
									&(cconf->kin_time));

						printf("Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s\n",
							cconf->exp_time, cconf->acc_time, cconf->kin_time);

						Sleep(1000);
					}

					// -----------------------------------
					//       query exposure time
					// -----------------------------------
					if (strncmp(cmd_in, "emgain ", 7) == 0)
					{
						sscanf(cmd_in, "%s %d", auxstr, &iChoice);
						AbortAcquisition(); // camera command

						error = SetEMCCDGain(iChoice);
						Sleep(1000);
					}

					// -----------------------------------
					//       query exposure time
					// -----------------------------------
					if (strncmp(cmd_in, "tint?", 5) == 0)
					{
						std::sprintf(msgout, "%.5f", cconf->exp_time);
						if (WriteFile(hPipeIn,
							msgout,
							std::strlen(msgout),   // = length of string + terminating '\0' !!!
							&dwWritten,
							NULL) == FALSE)
						{
							printf("ERROR: WriteFile(hPipeIn) failed!\n");
						}
						FlushFileBuffers(hPipeIn);
						//wfd = open(myfifout, O_WRONLY);
						//sprintf(msgout, "%.5f", cconf->exp_time);
						//write(wfd, msgout, LINESZ);
						//close(wfd);
						//wfd = 0;
					}
				}
				else
				{
					//DWORD dw = GetLastError(); 
					//printf("failed with error %d\n",dw);

					DisconnectNamedPipe(hPipeIn);
					CloseHandle(hPipeIn);
					hPipeIn = INVALID_HANDLE_VALUE;

					//if (dw == 109)
					//{
						hPipeIn = CreateNamedPipe(TEXT(myfifoin),
						PIPE_ACCESS_DUPLEX,
						PIPE_TYPE_MESSAGE | PIPE_READMODE_BYTE | PIPE_WAIT,   // FILE_FLAG_FIRST_PIPE_INSTANCE is not needed but forces CreateNamedPipe(..) to fail if the pipe already exists...
						1,
						1024 * 64,
						1024 * 64,
						PipeWait,
						NULL);

						//dw = GetLastError(); 
						//printf("failed with error %d\n",dw);
					//}

					Sleep(1000);
					break;
				}
			}
		}

		DisconnectNamedPipe(hPipeIn);
		DisconnectNamedPipe(hPipeOut);
	}
#endif

  printf("Shutting things down!\n");

  fflush(stdout);

  free(imsize);
  free(imarray);

  ShutDown();                // Shut down the iXon
  free(cconf);
  cconf = NULL;				// free cam_config structure
  return(EXIT_SUCCESS);
}

// --------------------------------------------------------------------------
//                              FUNCTION DEFINITIONS
// --------------------------------------------------------------------------
int CameraSelect (int iNumArgs, char* szArgList[]) {
  if (iNumArgs == 2) {
    
    at_32 lNumCameras;
    GetAvailableCameras(&lNumCameras);
    int iSelectedCamera = atoi(szArgList[1]);
    
    if (iSelectedCamera < lNumCameras && iSelectedCamera >= 0) {
      at_32 lCameraHandle;
      GetCameraHandle(iSelectedCamera, &lCameraHandle);
      SetCurrentCamera(lCameraHandle);
      return iSelectedCamera;
    }
    else
      return -1;
  }
  return(EXIT_SUCCESS);
}


int	StartThread_Acquisition(ThreadParams_Acquisition *param)
{
	param->hThreadCreated = CreateEvent(NULL,FALSE,FALSE,NULL);
	param->hTerminateThread = CreateEvent(NULL,FALSE,FALSE,NULL);
	param->hStopAcquisition = CreateEvent(NULL,FALSE,FALSE,NULL);
	param->hStartAcquisition = CreateEvent(NULL,FALSE,FALSE,NULL);
	param->hAcquisitionDone = CreateEvent(NULL,TRUE,TRUE,NULL);
	param->bEndThread = false;

	hThreadAcquisition = CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)Thread_Acquisition,(void*)param,THREAD_PRIORITY_ABOVE_NORMAL,&ThreadAcquisitionID);

	if (hThreadAcquisition == NULL)
	{
		// ERROR
		return -1;
	}

	int waitResult = WaitForSingleObject(param->hThreadCreated, THREAD_CREATION_TIMEOUT);

	switch(waitResult)
	{
	case WAIT_OBJECT_0:
		// THREAD CREATION SUCCESS
		printf("Thread creation SUCCESS!\n");
		break;
	case WAIT_FAILED:
		printf("Thread creation FAILED (WAIT_FAILED).\n");
		break;
	case WAIT_ABANDONED:
		printf("Thread creation FAILED (WAIT_ABANDONED).\n");
		break;
	case WAIT_TIMEOUT:
		printf("Thread creation FAILED (WAIT_TIMEOUT).\n");
		break;
	}

	return 0;
}

UINT32 Thread_Acquisition(void *pParam)
{
	ThreadParams_Acquisition *param = (ThreadParams_Acquisition*)pParam;
	
	HANDLE vEvent[3];
	vEvent[0] = param->hTerminateThread;
	vEvent[1] = param->hStopAcquisition;
	vEvent[2] = param->hStartAcquisition;

	int waitResult;

	SetEvent(param->hThreadCreated);

	while (param->bEndThread == false)
	{
		waitResult = WaitForMultipleObjects(3,&vEvent[0],FALSE,INFINITE);

		switch(waitResult)
		{
		case WAIT_OBJECT_0:
			// TERMINATE THREAD
			param->bEndThread = true;
			break;
		case WAIT_OBJECT_0+1:
			// STOP ACQUISITION
			//SetEvent(param->hAcquisitionStopped);
			break;
		case WAIT_OBJECT_0+2:
			// START ACQUISITION
			ResetEvent(param->hAcquisitionDone);
			acquire(param->pCamConfig);
			SetEvent(param->hAcquisitionDone);
			break;
		case WAIT_TIMEOUT:
			printf("Thread creation FAILED (WAIT_TIMEOUT).\n");
			break;
		}
	}

	return 0;
}