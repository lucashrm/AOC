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

#include <direct.h>

//#undef _TIMESPEC_DEFINED
#include <pthread.h>

#include "ImageStruct.h"
#include "ImageCreate.h"

#include <time.h>
#include <process.h>

#pragma comment(lib,"Advapi32.lib")

#include "AOC_server.h"

#include "thread_DMcmd_params.h"

// ANDOR TECHNOLOGY HEADER
#include <ATMCD32D.H>

// ALPAO HEADERS
#undef SUCCESS
#include <asdkDM.h>


#if defined ( WIN64 )
	#pragma comment(lib,"cfitsio_x64.lib")
	#pragma comment(lib,"libfftw-3.3_x64_dll.lib")
	#pragma comment(lib,"libcpuid_x64.lib")
	#pragma comment(lib,"mbgdevio_x64.lib")
	#pragma comment(lib,"mbgsvcio_x64.lib")
	#pragma comment(lib,"mbgutil_x64.lib")
	#pragma comment(lib,"gsl_x64.lib")
	#pragma comment(lib,"libeay32_x64.lib")
	#pragma comment(lib,"ssleay32_x64.lib")
	#if defined(_DEBUG)
		#pragma comment(linker, "/nodefaultlib:\"mfc100ud.lib\"")
		#pragma comment(linker, "/nodefaultlib:\"mfcs100ud.lib\"")
		#pragma comment(linker, "/nodefaultlib:\"libcmt.lib\"")
		#if defined(CUT_SECURE_SOCKET)
			#pragma comment(lib,"UTSecureLayer_x64_D.lib")
			#pragma comment(lib,"UltimateTCPIP_SS_x64_D.lib")
		#else
			#pragma comment(lib,"UltimateTCPIP_x64_D.lib")
		#endif
		#pragma comment(lib,"opencv_core231d_x64.lib")
		#pragma comment(lib,"opencv_imgproc231d_x64.lib")
		#pragma comment(lib,"opencv_objdetect231d_x64.lib")
		#pragma comment(lib,"opencv_gpu231d_x64.lib")
		#pragma comment(lib,"opencv_highgui231d_x64.lib")
		#pragma comment(lib,"CSMTP_x64_D.lib")
		#pragma comment(lib,"libcurl_openssl_libssh2_x64_D_dll.lib")
		#pragma comment(lib,"MCODE_Module_x64_D.lib")
		#pragma comment(lib,"MCODE_Data_x64_D.lib")
		#pragma comment(lib,"MCODE_Camera_x64_D.lib")
		#pragma comment(lib,"MCODE_Camera_Andor_x64_D.lib")
		#pragma comment(lib,"MCODE_LogManager_x64_D.lib")
		#pragma comment(lib,"MCODE_Communications_x64_D.lib")
		#pragma comment(lib,"MCODE_ConfigFileParser_x64_D.lib")
		#pragma comment(lib,"MCODE_Coordinates_x64_D.lib")
		#pragma comment(lib,"MCODE_Export_x64_D.lib")
		#pragma comment(lib,"MCODE_Export_Bitmap_x64_D.lib")
		#pragma comment(lib,"MCODE_Export_FITS_x64_D.lib")
		#pragma comment(lib,"MCODE_Ephem_x64_D.lib")
		#pragma comment(lib,"MCODE_GPS_x64_D.lib")
		#pragma comment(lib,"MCODE_GPS_Meinberg_x64_D.lib")
		#pragma comment(lib,"MCODE_Mail_x64_D.lib")
		#pragma comment(lib,"MCODE_MiscUtilities_x64_D.lib")
		#pragma comment(lib,"MCODE_SciLib_x64_D.lib")
		#pragma comment(lib,"MCODE_SciLib_OIS_x64_D.lib")
		#pragma comment(lib,"MCODE_SciLib_Proper_x64_D.lib")
		#pragma comment(lib,"MCODE_StdFunc_x64_D.lib")
		#pragma comment(lib,"MCODE_SystemUtilities_x64_D.lib")
		#pragma comment(lib,"MCODE_TimeSync_x64_D.lib")
		#pragma comment(lib,"MCODE_TimeSync_NTP_x64_D.lib")
	#else
		#pragma comment(linker, "/nodefaultlib:\"libcmt.lib\"")
		#pragma comment(lib,"libcurl_openssl_libssh2_x64_dll.lib")
		#if defined(CUT_SECURE_SOCKET)
			#pragma comment(lib,"UTSecureLayer_x64_D.lib")
			#pragma comment(lib,"UltimateTCPIP_SS_x64_D.lib")
		#else
			#pragma comment(lib,"UltimateTCPIP_x64_D.lib")
		#endif
		#pragma comment(lib,"opencv_core231_x64.lib")
		#pragma comment(lib,"opencv_imgproc231_x64.lib")
		#pragma comment(lib,"opencv_objdetect231_x64.lib")
		#pragma comment(lib,"opencv_gpu231_x64.lib")
		#pragma comment(lib,"opencv_highgui231_x64.lib")
		#pragma comment(lib,"CSMTP_x64.lib")
		#pragma comment(lib,"libcurl_openssl_libssh2_x64_dll.lib")
		#pragma comment(lib,"MCODE_Module_x64.lib")
		#pragma comment(lib,"MCODE_Data_x64.lib")
		#pragma comment(lib,"MCODE_Camera_x64.lib")
		#pragma comment(lib,"MCODE_Camera_Andor_x64.lib")
		#pragma comment(lib,"MCODE_LogManager_x64.lib")
		#pragma comment(lib,"MCODE_Communications_x64.lib")
		#pragma comment(lib,"MCODE_ConfigFileParser_x64.lib")
		#pragma comment(lib,"MCODE_Coordinates_x64.lib")
		#pragma comment(lib,"MCODE_Export_x64.lib")
		#pragma comment(lib,"MCODE_Export_Bitmap_x64.lib")
		#pragma comment(lib,"MCODE_Export_FITS_x64.lib")
		#pragma comment(lib,"MCODE_Ephem_x64.lib")
		#pragma comment(lib,"MCODE_GPS_x64.lib")
		#pragma comment(lib,"MCODE_GPS_Meinberg_x64.lib")
		#pragma comment(lib,"MCODE_Mail_x64.lib")
		#pragma comment(lib,"MCODE_MiscUtilities_x64.lib")
		#pragma comment(lib,"MCODE_SciLib_x64.lib")
		#pragma comment(lib,"MCODE_SciLib_OIS_x64.lib")
		#pragma comment(lib,"MCODE_SciLib_Proper_x64.lib")
		#pragma comment(lib,"MCODE_StdFunc_x64.lib")
		#pragma comment(lib,"MCODE_SystemUtilities_x64.lib")
		#pragma comment(lib,"MCODE_TimeSync_x64.lib")
		#pragma comment(lib,"MCODE_TimeSync_NTP_x64.lib")
	#endif
#elif defined ( WIN32 )
	#pragma comment(lib,"cfitsio.lib")

	#pragma comment(lib,"libeay32.lib")
	#pragma comment(lib,"ssleay32.lib")


	#pragma comment(lib,"ATMCD32M.LIB")
	#pragma comment(lib,"libfftw-3.3_x86_dll.lib")
	#pragma comment(lib,"mbgdevio.lib")
	#pragma comment(lib,"mbgsvcio.lib")
	#pragma comment(lib,"mbgutil.lib")
	#pragma comment(lib,"gsl.lib")
	#pragma comment(lib,"ftd2xx_static_x86.lib")
	#if defined(_DEBUG)
		#if defined(CUT_SECURE_SOCKET)
			#pragma comment(lib,"UTSecureLayer_D.lib")
			#pragma comment(lib,"UltimateTCPIP_SS_D.lib")
		#else
			#pragma comment(lib,"UltimateTCPIP_D.lib")
		#endif
		#pragma comment(lib,"CSMTP_D.lib")
		#pragma comment(lib,"TIS_UDSHL10d.lib")
		#pragma comment(lib,"libcurl_openssl_libssh2_D_dll.lib")
		//#pragma comment(linker, "/nodefaultlib:\"msvcrtd.lib\"")
		#pragma comment(linker, "/nodefaultlib:\"libcmt.lib\"")
		//#pragma comment(linker, "/nodefaultlib:\"libcmtd.lib\"")
		#pragma comment(lib,"CGridListCtrlEx_D.lib")
		#pragma comment(lib,"MCODE_Module_D.lib")
		#pragma comment(lib,"MCODE_Data_D.lib")
		#pragma comment(lib,"MCODE_Camera_D.lib")
		#if defined ( _MCODE_USE_CAMERA_APOGEE )
			#pragma comment(lib,"MCODE_Camera_Apogee_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_SBIG )
			#pragma comment(lib,"MCODE_Camera_SBIG_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_QSI )
			#pragma comment(lib,"MCODE_Camera_QSI_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_ANDOR )
			#pragma comment(lib,"MCODE_Camera_Andor_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_FLI )
			#pragma comment(lib,"MCODE_Camera_FLI_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_FLIPRO )
			#pragma comment(lib,"libflipro.lib")
			#pragma comment(lib,"MCODE_Camera_FLIPRO_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_IC )
			#pragma comment(lib,"MCODE_Camera_IC_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_INOVA )
			#pragma comment(lib,"libinovasdkd.lib")
			#pragma comment(lib,"MCODE_Camera_iNova_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_RAPTORPHOTONICS )
			#pragma comment(lib,"MCODE_Camera_RaptorPhotonics_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_BASLER )
			#pragma comment(lib,"MCODE_Camera_Basler_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_QHY )
			#pragma comment(lib,"MCODE_Camera_QHY_D.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_ATIK )
			#pragma comment(lib,"MCODE_Camera_Atik_D.lib")
		#endif
		#pragma comment(lib,"MCODE_LogManager_D.lib")
		#pragma comment(lib,"MCODE_DeformableMirror_D.lib")
		#pragma comment(lib,"MCODE_DeformableMirror_BMC_D.lib")
		#pragma comment(lib,"MCODE_DLL_LogManagerGUI_v02_D.lib")
		#pragma comment(lib,"MCODE_RichEdit_D.lib")
		#pragma comment(lib,"MCODE_CmdLineEditCtrl_D.lib")
		#pragma comment(lib,"MCODE_Communications_D.lib")
		#pragma comment(lib,"MCODE_GraphxDisplay_D.lib")
		#pragma comment(lib,"MCODE_Mail_D.lib")
		#pragma comment(lib,"MCODE_MiscUtilities_D.lib")
		#pragma comment(lib,"MCODE_Mount_D.lib")
		#pragma comment(lib,"MCODE_Mount_AP_D.lib")
		#pragma comment(lib,"MCODE_Mount_C2PUTCS_D.lib")
		#pragma comment(lib,"MCODE_Export_D.lib")
		#pragma comment(lib,"MCODE_Export_Bitmap_D.lib")
		#pragma comment(lib,"MCODE_Export_FITS_D.lib")
		#pragma comment(lib,"MCODE_TimeSync_D.lib")
		#pragma comment(lib,"MCODE_TimeSync_NTP_D.lib")
		#pragma comment(lib,"MCODE_GPS_D.lib")
		#pragma comment(lib,"MCODE_GPS_Meinberg_D.lib")
		#pragma comment(lib,"MCODE_Coordinates_D.lib")
		#pragma comment(lib,"MCODE_CondSolver_D.lib")
		#pragma comment(lib,"MCODE_Dome_D.lib")
		#pragma comment(lib,"MCODE_Dome_AshDome_D.lib")
		#pragma comment(lib,"MCODE_FilterWheel_D.lib")
		#pragma comment(lib,"MCODE_FilterWheel_Thorlabs_D.lib")
		#pragma comment(lib,"MCODE_GraphxLabel_D.lib")
		#pragma comment(lib,"MCODE_GraphScheduler_D.lib")
		#pragma comment(lib,"MCODE_ConfigFileParser_D.lib")
		#pragma comment(lib,"MCODE_ModuleOptionPanel_D.lib")
		//#pragma comment(lib,"MCODE_MotController_D.lib")
		//#pragma comment(lib,"MCODE_MotController_Newport_D.lib")
		#pragma comment(lib,"MCODE_NamedVariable_D.lib")
		#pragma comment(lib,"MCODE_SystemUtilities_D.lib")
		#pragma comment(lib,"MCODE_Graph_D.lib")
		#pragma comment(lib,"MCODE_Ephem_D.lib")
		#pragma comment(lib,"MCODE_Focuser_D.lib")
		#pragma comment(lib,"MCODE_Focuser_Micos_D.lib")
		#pragma comment(lib,"MCODE_Focuser_MicosHydra_D.lib")
		#pragma comment(lib,"MCODE_Focuser_Optec_D.lib")
		#pragma comment(lib,"MCODE_RemotePwr_D.lib")
		#pragma comment(lib,"MCODE_RemotePwr_BBPS_D.lib")
		#pragma comment(lib,"MCODE_SchedCond_D.lib")
		#pragma comment(lib,"MCODE_SchedCond_SunHeight_D.lib")
		#pragma comment(lib,"MCODE_SchedCond_TimeBounds_D.lib")
		#pragma comment(lib,"MCODE_SchedCond_Periodic_D.lib")
		#pragma comment(lib,"MCODE_SciLib_D.lib")
		#pragma comment(lib,"MCODE_SciLib_OIS_D.lib")
		#pragma comment(lib,"MCODE_SciLib_Proper_D.lib")
		#pragma comment(lib,"MCODE_Serial_D.lib")
		#pragma comment(lib,"MCODE_SpinCounter_D.lib")
		#pragma comment(lib,"MCODE_SpinCounter_Raspberry_D.lib")
		#pragma comment(lib,"MCODE_StarCatalog_D.lib")
		#pragma comment(lib,"MCODE_StdFunc_D.lib")
	#else
		#pragma comment(linker, "/nodefaultlib:\"libcmt.lib\"")
		#pragma comment(lib,"CGridListCtrlEx.lib")
		#pragma comment(lib,"TIS_UDSHL10.lib")
		#pragma comment(lib,"libcurl_openssl_libssh2_dll.lib")
		#if defined(CUT_SECURE_SOCKET)
			#pragma comment(lib,"UTSecureLayer.lib")
			#pragma comment(lib,"UltimateTCPIP_SS.lib")
		#else
			#pragma comment(lib,"UltimateTCPIP.lib")
		#endif
		#pragma comment(lib,"CSMTP.lib")
		#pragma comment(lib,"MCODE_Camera.lib")
		#if defined ( _MCODE_USE_CAMERA_APOGEE )
			#pragma comment(lib,"MCODE_Camera_Apogee.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_SBIG )
			#pragma comment(lib,"MCODE_Camera_SBIG.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_QSI )
			#pragma comment(lib,"MCODE_Camera_QSI.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_ANDOR )
			#pragma comment(lib,"MCODE_Camera_Andor.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_FLI )
			#pragma comment(lib,"libflipro.lib")
			#pragma comment(lib,"MCODE_Camera_FLI.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_FLIPRO )
			#pragma comment(lib,"MCODE_Camera_FLIPRO.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_IC )
			#pragma comment(lib,"MCODE_Camera_IC.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_INOVA )
			#pragma comment(lib,"MCODE_Camera_iNova.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_RAPTORPHOTONICS )
			#pragma comment(lib,"MCODE_Camera_RaptorPhotonics.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_BASLER )
			#pragma comment(lib,"MCODE_Camera_Basler.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_QHY )
			#pragma comment(lib,"MCODE_Camera_QHY.lib")
		#endif
		#if defined ( _MCODE_USE_CAMERA_ATIK )
			#pragma comment(lib,"MCODE_Camera_Atik.lib")
		#endif
		#pragma comment(lib,"MCODE_CmdLineEditCtrl.lib")
		#pragma comment(lib,"MCODE_Communications.lib")
		#pragma comment(lib,"MCODE_Data.lib")
		#pragma comment(lib,"MCODE_DeformableMirror.lib")
		#pragma comment(lib,"MCODE_DeformableMirror_BMC.lib")
		#pragma comment(lib,"MCODE_DLL_LogManagerGUI_v02.lib")
		#pragma comment(lib,"MCODE_Ephem.lib")
		#pragma comment(lib,"MCODE_LogManager.lib")
		#pragma comment(lib,"MCODE_Module.lib")
		#pragma comment(lib,"MCODE_RichEdit.lib")
		#pragma comment(lib,"MCODE_GraphxDisplay.lib")
		#pragma comment(lib,"MCODE_Mail.lib")
		#pragma comment(lib,"MCODE_MiscUtilities.lib")
		#pragma comment(lib,"MCODE_Mount.lib")
		#pragma comment(lib,"MCODE_Mount_AP.lib")
		#pragma comment(lib,"MCODE_Mount_C2PUTCS.lib")
		#pragma comment(lib,"MCODE_Export.lib")
		#pragma comment(lib,"MCODE_Export_Bitmap.lib")
		#pragma comment(lib,"MCODE_Export_FITS.lib")
		#pragma comment(lib,"MCODE_TimeSync.lib")
		#pragma comment(lib,"MCODE_TimeSync_NTP.lib")
		#pragma comment(lib,"MCODE_GPS.lib")
		#pragma comment(lib,"MCODE_GPS_Meinberg.lib")
		#pragma comment(lib,"MCODE_Coordinates.lib")
		#pragma comment(lib,"MCODE_CondSolver.lib")
		#pragma comment(lib,"MCODE_Dome.lib")
		#pragma comment(lib,"MCODE_Dome_AshDome.lib")
		#pragma comment(lib,"MCODE_GraphxLabel.lib")
		#pragma comment(lib,"MCODE_GraphScheduler.lib")
		#pragma comment(lib,"MCODE_FilterWheel.lib")
		#pragma comment(lib,"MCODE_FilterWheel_Thorlabs.lib")
		#pragma comment(lib,"MCODE_ConfigFileParser.lib")
		#pragma comment(lib,"MCODE_ModuleOptionPanel.lib")
		//#pragma comment(lib,"MCODE_MotController.lib")
		//#pragma comment(lib,"MCODE_MotController_Newport.lib")
		#pragma comment(lib,"MCODE_NamedVariable.lib")
		#pragma comment(lib,"MCODE_SystemUtilities.lib")
		#pragma comment(lib,"MCODE_Graph.lib")
		#pragma comment(lib,"MCODE_Focuser.lib")
		#pragma comment(lib,"MCODE_Focuser_Micos.lib")
		#pragma comment(lib,"MCODE_Focuser_MicosHydra.lib")
		#pragma comment(lib,"MCODE_Focuser_Optec.lib")
		#pragma comment(lib,"MCODE_RemotePwr.lib")
		#pragma comment(lib,"MCODE_RemotePwr_BBPS.lib")
		#pragma comment(lib,"MCODE_SchedCond.lib")
		#pragma comment(lib,"MCODE_SchedCond_SunHeight.lib")
		#pragma comment(lib,"MCODE_SchedCond_TimeBounds.lib")
		#pragma comment(lib,"MCODE_SchedCond_Periodic.lib")
		#pragma comment(lib,"MCODE_SciLib.lib")
		#pragma comment(lib,"MCODE_SciLib_OIS.lib")
		#pragma comment(lib,"MCODE_SciLib_Proper.lib")
		#pragma comment(lib,"MCODE_Serial.lib")
		#pragma comment(lib,"MCODE_SpinCounter.lib")
		#pragma comment(lib,"MCODE_SpinCounter_Raspberry.lib")
		#pragma comment(lib,"MCODE_StarCatalog.lib")
		#pragma comment(lib,"MCODE_StdFunc.lib")
	#endif
#endif

// PLATFORMLESS LIBRARIES 
#pragma comment(lib,"ASDK.lib")

using namespace std;
using namespace MultiCODE;
using namespace acs;

#define MAX_BUFFER_HISTORY_LENGTH		int(500)


#define LINESZ 256
char myfifoin[LINESZ] = "\\\\.\\pipe\\ixon_fifo_in";  // incoming pipe
char myfifout[LINESZ] = "\\\\.\\pipe\\ixon_fifo_out"; // out-going pipe

static string REVISION = "039";
MCODE_LogManager *Log;

// MultiCODE LOG
string				LogMsg;

// DEFORMABLE MIRROR INSTANCE
#define	DM_ID		"BOL115"
//acs::DM				*pDM = NULL;
// DM INSTANCE CREATION
acs::DM				dm(DM_ID);

bool				bSimu_iXon;
bool				bSimu_DM;
int					nAct;
int					DM_nAct;

// IXON CAMERA THREAD PARAMETERS
ThreadParams_iXon	TP_iXon;
DWORD				Thread_iXon_ID;
HANDLE				hThread_iXon;

// SIMULATION IXON THREAD PARAMETERS
ThreadParams_Simu_iXon	TP_Simu_iXon;
DWORD					Thread_Simu_iXon_ID;
HANDLE					hThread_Simu_iXon;


// ALPAO DM THREAD PARAMETERS
ThreadParams_DMcmd	TP_DMcmd;
DWORD				Thread_DMcmd_ID;
HANDLE				hThread_DMcmd;

// SHACK-HARTMANN CALCULUS THREAD PARAMETERS
ThreadParams_calcSH	TP_calcSH;
DWORD 				Thread_calcSH_ID;
HANDLE				hThread_calcSH;

// EVENT LISTENER THREAD PARAMETERS
ThreadParams_EvtWatch	TP_EvtWatch;
DWORD 				Thread_EvtWatch_ID;
HANDLE				hThread_EvtWatch;

// DATA LOGGER THREAD PARAMETERS
ThreadParams_DataLogger	TP_DataLogger;
DWORD 					Thread_DataLogger_ID;
HANDLE					hThread_DataLogger;

// DATA SAVER THREAD PARAMETERS (SEE DATA LOGGER IMPLEMENTATION FILE)
ThreadParams_DataSaver	TP_DataSaver;
DWORD 					Thread_DataSaver_ID;
HANDLE					hThread_DataSaver;


IMAGE				*SHM_DM_state;			// SHARED MEMORY DM STATE
IMAGE				*SHM_imarray;			// SHARED MEMORY INSTANT IXON IMAGE
IMAGE				*SHM_Slopes;			// SHARED MEMORY N FRAMES HISTORY DATA (DISPLAY RATE)
IMAGE				*SHM_RefCent;			// SHARED MEMORY REFERENCE CENTROIDS DATA
IMAGE				*SHM_DM_slope_mask;		// SHARED MEMORY PUPIL MASK DATA
IMAGE				*SHM_ModesInfo;			// SHARED MEMORY MODE INFORMATION
IMAGE				*SHM_Timestamps;        // SHARED MEMORY TIMESTAMPS
vector<IMAGE*>		vSHM_DM_BUFFERS;		// SHARED MEMORY CONTAINING DM MAPS SEPARATE CHANNELS

vector<acs::Scalar>	vDM_cmd;

// INTERNAL DM COMMAND MAPS (nAct SIZE)
int						N_DM_COMMAND_MAPS = 8;					// NUMBER OF CORRECTION BUFFERS (NOT INCLUDING THE DM COLLAPSED BUFFER => SHM_DM_state)

vector<vector<float>>	vDM_cmd_buffers;		// DM CORRECTION BUFFERS (N_DM_COMMAND_MAPS x nAct)
vector<float>			vDM_static_cmd;			// DM STATIC COMMAND MAP (MAP CORRESPONDING TO NON-TIME-CRITICAL MAPS)
vector<size_t>			vCollapseStaticInd;		// DM STATIC MAPS INDICES
vector<size_t>			vCollapseDynamicInd;	// DM DYNAMIC (TIME-CRITICAL) INDICES = NOT(vCollapseStaticInd)
vector<float>			vDM_cur_cmd;			// DM INSTANTANEOUS COMMAND MAP = vDM_cmd_buffers[vCollapseIndices] + vDM_static_cmd
vector<long long>		vTimeStamp;				// TIME VECTOR RECORDED JUST BEFORE SETTING THE DM INTO POSITION

vector<float>			vtmpLO;
vector<float>			vtmpHO;
vector<double>			vTT_drift;
vector<double>			vHO_drift;

vector<float>			vCmdMatrix;

float					leakage_gain;

vector<int>				vDMSquare2DMmap;

int						localCount;

mcINT16 concat_quoted_string_items(string& output, vector<string> vItem, int imin = -1, int imax = -1, bool bKeepQuotes = false)
{
	if (imin == -1)
	{
		imin = 0;
	}

	if (imax == -1)
	{
		imax = int(vItem.size()) - 1;
	}

	std::string strtmp;
	size_t posa;
	size_t posb;

	bool bFound = false;

	while(imin <= imax)
	{
		if ((posa = vItem[imin].find("\"")) != string::npos)
		{
			if ((posb = vItem[imin].find("\"",posa+1)) != string::npos)
			{
				strtmp = vItem[imin].substr(posa, posb - posa + 1);
				bFound = true;
				break;
			}
			else
			{
				strtmp = vItem[imin].substr(posa);

				for (size_t k=imin+1;k<=imax;++k)
				{
					posa = vItem[k].find("\"");
					if (posa == string::npos)
					{
						strtmp += vItem[k];
						continue;
					}
					else
					{
						strtmp += vItem[k].substr(0, posa+1);
						bFound = true;
						break;
					}
				}

				if (bFound == true)
				{
					break;
				}
			}
		}
		else
		{
			if ((imin < vItem.size()) && (imin <= imax))
			{
				imin += 1;
			}
		}
	}

	if (bFound == false)
	{
		return WARNING;
	}

	if (bKeepQuotes == true)
	{
		output = strtmp;
	}
	else
	{
		output = strtmp.substr(1, strtmp.size() - 2);
	}
	return SUCCESS;
}

mcINT16 ApplyRefCentroids()
{
	// SET THE X AND Y REFERENCE CENTROID VALUES									
	std::copy(
		(float*)SHM_RefCent[0].array.SI32,
		(float*)SHM_RefCent[0].array.SI32 + TP_calcSH.nx*TP_calcSH.ny,
		&TP_calcSH.vSH_ref_x[0]);

	std::copy(
		(float*)SHM_RefCent[0].array.SI32 + TP_calcSH.nx*TP_calcSH.ny,
		(float*)SHM_RefCent[0].array.SI32 + 2*TP_calcSH.nx*TP_calcSH.ny,
		&TP_calcSH.vSH_ref_y[0]);

	return SUCCESS;
}

mcINT16 shutdown_dm()
{
	if (dm.Check() == true)
	{
		vector<acs::Scalar> zero_map(DM_nAct, 0.);
		dm.Send(zero_map.data());
		dm.Stop();
	}

	return SUCCESS;
}

void main_loop()
{
	char *cmd_in = (char*) malloc(LINESZ * sizeof(char));
	char *buf = (char*) malloc(LINESZ * sizeof(char));
	char *auxstr = (char*) malloc(LINESZ * sizeof(char));
	char *msgout = (char*) malloc(LINESZ * sizeof(char));

	MCODE_Export_FITS FITS;

	// ------------------------------------------------
	// create pipes for interaction with other programs
	// ------------------------------------------------
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

	// ------------------------------------------
	//          main server loop
	// ------------------------------------------
	int rfd = 0; // file descriptor for input fifo
	int wfd = 0; // file descriptor for output fifo
	int foo = -1;
	bool bQuit = false;
	//int iChoice;
	//float fChoice;
	unsigned int state;
	mcUINT32 error;

	//ThreadParams_Acquisition TP_ThreadAcquisition;
	//TP_ThreadAcquisition.pCamConfig = cconf;
	//
	//if (StartThread_Acquisition(&TP_ThreadAcquisition) == 0)
	//{
	//	printf("Thread_Acquisition started SUCCESSFULLY.\n");
	//	quit = false;
	//}
	//else
	//{
	//	printf("Thread_Acquisition FAILED to start.\n");
	//	quit = true;
	//}

	HANDLE hWaitEvent = CreateEvent(NULL,FALSE,FALSE,NULL);

	cout << string_format("\niXon server ready.\n");
	//std::fflush(stdout);

	DWORD dwRead;

	//hPipeIn = INVALID_HANDLE_VALUE;

	if (hPipeIn != INVALID_HANDLE_VALUE)
	{

		while((bQuit == false) || (ConnectNamedPipe(hPipeIn, NULL) != FALSE))   // wait for someone to connect to the pipe
		{
			int waitResult = WaitForSingleObject(hWaitEvent,1000);
			//DWORD dw = GetLastError(); 
			//printf("failed with error %d\n",dw);
			while (bQuit == false)
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

					string myCmd = cmd_in;

					cout << string_format("Received command: \"%s\"\n", myCmd.data());

					trim(myCmd," ");
					vector<string> vItem = split(myCmd,' ');

					if (vItem.size() < 2)
					{
						cout << string_format("Incorrect number of arguments (syntaxe: <module> <command> [<parameters>])\n", myCmd.data());
						continue;
					}

					string module = vItem[0];
					string cmdStr = vItem[1];

					while(true)
					{
						if (compare_no_case(module,"ixon") == 0)
						{
							/* do something with data in buffer */
							// ===============================================
							//       processing the command from pipe
							// ===============================================

							// -----------------------------------
							//       continuous acquisition
							// -----------------------------------
							if (compare_no_case(cmdStr, "stream") == 0)
							{
								if (TP_iXon.bStreamON == false)
								{
									//// SET THE X AND Y REFERENCE CENTROID VALUES
									//std::copy(
									//	(float*)SHM_RefCent[0].array.SI32,
									//	(float*)SHM_RefCent[0].array.SI32 + TP_calcSH.nx*TP_calcSH.ny,
									//	&TP_calcSH.vSH_ref_x[0]);

									//std::copy((float*)SHM_RefCent[0].array.SI32 + TP_calcSH.nx*TP_calcSH.ny,
									//	(float*)SHM_RefCent[0].array.SI32 + 2*TP_calcSH.nx*TP_calcSH.ny,
									//	&TP_calcSH.vSH_ref_y[0]);

									// SET THE X AND Y REFERENCE CENTROID VALUES
									ApplyRefCentroids();

									// INITIALIZE DISPLAY BUFFERS
									TP_calcSH.display_image_count = TP_iXon.display_image_count;
									//TP_calcSH.vFlatSHData.resize(TP_calcSH.display_image_count * TP_calcSH.nx * TP_calcSH.ny * TP_calcSH.nDataSets);
									TP_calcSH.vModesDataHist.resize(TP_calcSH.display_image_count*TP_calcSH.MaxModes);

									vector<float> vtmp(TP_calcSH.display_image_count*TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets, 0.f);
									vector<double> vtimestamptmp(TP_calcSH.display_image_count * 50, 0);
									vector<float> vdmtmp(TP_calcSH.display_image_count * 97, 0);
									TP_calcSH.vCircBuf_FlatSHData.resize(TP_calcSH.nCircBuffers, vtmp);
									TP_calcSH.vCircBuf_Timestamps.resize(TP_calcSH.nCircBuffers, vtimestamptmp);
									TP_calcSH.vCircBuf_DM_cmd.resize(TP_calcSH.nCircBuffers, vdmtmp);

									TP_calcSH.circBufCounter = 0;
									TP_calcSH.local_counter = 0;

									TP_iXon.bAcqON = true;
									TP_iXon.bStreamON = true;

									SetEvent(TP_Simu_iXon.hStartLiveStream);
									SetEvent(TP_iXon.hStartAcquisition);
								}

								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//      abort current acquisition
							// -----------------------------------
							else if (compare_no_case(cmdStr,"stop") == 0)
							{
								if (TP_iXon.bAcqON)
								{
									//cconf->babort = true;
									SetEvent(TP_iXon.hStopAcquisition);
									SetEvent(TP_Simu_iXon.hStopLiveStream);
									if (WaitForSingleObject(TP_iXon.hAcquisitionDone,5000) == WAIT_TIMEOUT)
									{
										cout << string_format("TIMEOUT while waiting for acquisition routine to terminate\n");
									}

									//::AbortAcquisition(); // camera command
									TP_iXon.bAcqON = false;
									TP_iXon.bStreamON = false;
								}
								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//  update the subaperture thresholds
							// -----------------------------------
							else if (compare_no_case(cmdStr,"set_thresh") == 0 && !bSimu_iXon)
							{
								if (vItem.size() < 4)
								{
									cout << string_format("ERROR: 'ixon set_thresh \"<strFile>\" \"<strFile>\"' Incorrect number of parameters (2 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								string strFileThresh;
								string strFileEnable;
								concat_quoted_string_items(strFileThresh, vItem, 2, 2);
								concat_quoted_string_items(strFileEnable, vItem, 3, 3);

								vector<int> vSize;
								vector<float> vSubApThresh;
								if (FITS.Load(strFileThresh, vSubApThresh, vSize) != SUCCESS)
								{
									FITS.CloseFITS();
									cout << string_format("ERROR: failed to load subaperture threshold FITS file \"%s\".\n", strFileThresh.data());
									break;
								}
								if (vSize[0] != TP_calcSH.nx*TP_calcSH.ny)
								{
									cout << string_format("ERROR: 'ixon set_thresh': subaperture thresholds array must contain %d float values. Nothing done.\n", TP_calcSH.nAct);
									break;
								}
								std::copy(vSubApThresh.begin(),vSubApThresh.end(),TP_calcSH.vThresh.begin());


								vector<int> vSubApEnable;
								if (FITS.Load(strFileEnable, vSubApEnable, vSize) != SUCCESS)
								{
									FITS.CloseFITS();
									cout << string_format("ERROR: failed to load subaperture enabled FITS file \"%s\".\n", strFileEnable.data());
									break;
								}
								if (vSize[0] != TP_calcSH.nx*TP_calcSH.ny)
								{
									cout << string_format("ERROR: 'ixon set_thresh': enabled subaperture array must contain %d float values. Nothing done.\n", TP_calcSH.nAct);
									break;
								}
								std::copy(vSubApEnable.begin(),vSubApEnable.end(),TP_calcSH.vEnabled.begin());

								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//  sets the camera transfer mode
							// -----------------------------------
							else if (compare_no_case(cmdStr,"set_transfer_mode") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("ERROR: 'ixon set_transfer_mode <int>' Incorrect number of parameters (1 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								int mode = std::atoi(vItem[2].data());
								
								if ((mode < 0) || (mode > 1))
								{
									cout << string_format("ERROR: 'ixon set_transfer_mode <int>' incorrect mode range (0=Normal or 1=FT). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								TP_iXon.FTmode = mode;
								::SetFrameTransferMode(TP_iXon.FTmode);
								AOC_SetHWConf(&TP_iXon);

								cout << string_format("iXon frame transfer mode set to %d (0=Normal or 1=FT).\n", mode);
								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//         close the program
							// -----------------------------------
							//if (strncmp(cmd_in, "quit", 4) == 0)
							else if (compare_no_case(cmdStr,"quit") == 0)
							{
								if (TP_iXon.bAcqON)
								{
									//cconf->babort = true;
									SetEvent(TP_iXon.hStopAcquisition);
									if (WaitForSingleObject(TP_iXon.hAcquisitionDone,5000) == WAIT_TIMEOUT)
									{
										cout << string_format("TIMEOUT while waiting for acquisition routine to terminate\n");
									}
									//AbortAcquisition(); // camera command
									cout << string_format("\nQuitting program\n");
									//Sleep(500);
								}

								if (hThread_DMcmd != INVALID_HANDLE_VALUE)
								{
									TerminateThread_DMcmd();
								}

								bQuit = true;
								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//       set temperature goal
							// -----------------------------------
							//if (strncmp(cmd_in, "tempg ", 6) == 0)
							else if (compare_no_case(cmdStr,"tempg") == 0)
							{
								//sscanf(cmd_in, "%s %d", auxstr, &iChoice);
								int temp = std::atoi(vItem[2].data());
								if ((temp < -60) || (temp > 20))
								{
									cout << string_format("\nTemperature set point (%d) not in range (-60<T<+20). Nothing done.\n",temp);
								}
								else
								{
									SetTemperature(temp); // camera command
									cout << string_format("\nTempertaure set point set to %d Cdeg.\n",temp);
								}
								//Sleep(200);
								break;		// TO EXIT THE while(true) STATEMENT
							}


							// -----------------------------------
							//       set temperature goal
							// -----------------------------------
							//if (strncmp(cmd_in, "cooling ", 8) == 0)
							else if (compare_no_case(cmdStr,"cooling") == 0)
							{
								//sscanf(cmd_in, "%s %d", auxstr, &iChoice);
								bool bOn = std::atoi(vItem[2].data()) == 0 ? false : true;
								if (bOn == false)
								{
									state = ::CoolerOFF();
									if (state != DRV_SUCCESS)
									{
										cout << string_format("\n'::CoolerOFF' FAILED to execute (error code=%d).\n",state);
									}
									else
									{
										cout << string_format("\nCooling is now OFF.\n");
									}
								}
								else if (bOn == true)
								{
									state = ::CoolerON();
									if (state != DRV_SUCCESS)
									{
										cout << string_format("\n'::CoolerON' FAILED to execute (error code=%d).\n",state);
									}
									else
									{
										cout << string_format("\nCooling is now ON.\n");
									}
								}
								//Sleep(200);
								break;		// TO EXIT THE while(true) STATEMENT
							}


							//if (strncmp(cmd_in, "fTemp?", 6) == 0)
							else if (compare_no_case(cmdStr,"fTemp?") == 0)
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
									cout << string_format("ERROR: WriteFile(hPipeIn) failed!\n");
								}
								FlushFileBuffers(hPipeIn);
	
								//Sleep(200);
								break;		// TO EXIT THE while(true) STATEMENT
							}

							//// -----------------------------------
							////       take a single image
							//// -----------------------------------
							//if (strncmp(cmd_in, "take ", 5) == 0)
							//{
							//	foo = -1;
							//	sscanf(cmd_in, "%s %d", auxstr, &foo);
							//	cconf->nleft = (foo < 0) ? 1: foo; 
							//	cconf->bacqon = true;
							//	pthread_create(&tid_acqr, NULL, acquire, cconf);
							//}


							// -----------------------------------
							//       update exposure time
							// -----------------------------------
							//if (strncmp(cmd_in, "tint ", 5) == 0)
							else if (compare_no_case(cmdStr,"tint") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								//sscanf(cmd_in, "%s %f", auxstr, &fChoice);
								//printf("Exposure time: %f\n", fChoice);
								float exptime = float(std::atof(vItem[2].data()));

								if (TP_iXon.bAcqON == true)
								{
									SetEvent(TP_iXon.hStopAcquisition);
									if (WaitForSingleObject(TP_iXon.hAcquisitionDone,5000) == WAIT_TIMEOUT)
									{
										cout << string_format("TIMEOUT while waiting for acquisition routine to terminate\n");
									}
								}
								
								TP_iXon.expTime_s = exptime;
								TP_iXon.kinTime_s = TP_iXon.expTime_s;
								SHM_imarray->kw[0].value.numl = int(exptime * 1000000.f + 0.5f);

								//SetFrameTransferMode(1);					// Frame transfer mode is OFF
								//SetReadMode(4);							// Set Read Mode to --Image--
								//SetAcquisitionMode(5);					// Set Acq. mode to --Run till abort--
								//SetExposureTime(TP_iXon.expTime_s);
								//SetKineticCycleTime(0.f);
								if (TP_iXon.FTmode == 0)
								{
									SetExposureTime(TP_iXon.expTime_s);
									SetKineticCycleTime(0.f);
								}
								else if (TP_iXon.FTmode == 1)
								{
									SetExposureTime(TP_iXon.expTime_s);
									SetKineticCycleTime(TP_iXon.kinTime_s);
								}

								AOC_SetHWConf(&TP_iXon);
								TP_calcSH.display_image_count = TP_iXon.display_image_count;
								//TP_calcSH.vFlatSHData.resize(TP_calcSH.display_image_count*TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets);
								TP_calcSH.vModesDataHist.resize(TP_calcSH.display_image_count*TP_calcSH.MaxModes);
								
								vector<float> vtmp(TP_calcSH.display_image_count*TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets, 0.f);
								TP_calcSH.vCircBuf_FlatSHData.resize(TP_calcSH.nCircBuffers, vtmp);

								//Sleep(200);
								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//       set EM CCD gain
							// -----------------------------------
							//if (strncmp(cmd_in, "emgain ", 7) == 0)
							else if (compare_no_case(cmdStr,"emgain") == 0)
							{
								//sscanf(cmd_in, "%s %d", auxstr, &iChoice);
								int gain = std::atoi(vItem[2].data());
								SetEvent(TP_iXon.hStopAcquisition);
								if (WaitForSingleObject(TP_iXon.hAcquisitionDone,5000) == WAIT_TIMEOUT)
								{
									cout << string_format("TIMEOUT while waiting for acquisition routine to terminate\n");
								}

								error = SetEMCCDGain(gain);
								Sleep(1000);
								break;		// TO EXIT THE while(true) STATEMENT
							}

							// -----------------------------------
							//       query exposure time
							// -----------------------------------
							//if (strncmp(cmd_in, "tint?", 5) == 0)
							else if (compare_no_case(cmdStr,"tint?") == 0)
							{
								std::sprintf(msgout, "%.6f", TP_iXon.expTime_s);
								if (WriteFile(hPipeIn,
									msgout,
									std::strlen(msgout),   // = length of string + terminating '\0' !!!
									&dwWritten,
									NULL) == FALSE)
								{
									cout << string_format("ERROR: WriteFile(hPipeIn) failed!\n");
								}
								FlushFileBuffers(hPipeIn);
								break;		// TO EXIT THE while(true) STATEMENT
							}

							break;		// TO EXIT THE while(true) STATEMENT
						}
						else if (compare_no_case(module,"alpao") == 0 && !bSimu_DM)
						{
							if (compare_no_case(cmdStr,"enable_dm") == 0 && bSimu_DM)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("'alpao enable_dm <0/1>' incorrect number of arguments. Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								TP_DMcmd.bEnableDM = std::atoi(vItem[2].data()) == 0 ? false : true;

								cout << string_format("DM realtime update %s.\n", TP_DMcmd.bEnableDM == true ? "ON" : "OFF");
							}
							else if (compare_no_case(cmdStr,"send_cmdmap_to_dm") == 0)
							{
								if (vItem.size() != 2)
								{
									cout << string_format("'alpao send_cmdmap_to_dm' incorrect number of arguments. Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								if (dm.Check() != true)
								{
									cout << string_format("DM pointer (pDM) is NULL! Nothing done.\n");
									break;
								}


								//vector<float>::iterator itf = vDM_cur_cmd.begin();
								//vector<double>::iterator itd = vDM_cmd.begin();
								//for (size_t k = 0; k < vDM_cmd.size(); ++k, ++itf, ++itd)
								//{
								//	*itd = double(*itf);
								//}
								//
								////dm.Reset();
								//dm.Send(vDM_cmd.data());

								DM_set_to_position();

								cout << string_format("DM position updated with combined command map.\n");
							}
							else if (compare_no_case(cmdStr,"update_static_maps") == 0)
							{
								if (vItem.size() != 2)
								{
									cout << string_format("'alpao update_combined' incorrect number of parameters. Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								//collapse_DM_buffers();
								update_DM_static_buffer();

								cout << string_format("'alpao update_static_maps' command issued.\n");
							}
							else if (compare_no_case(cmdStr,"pulse") == 0)
							{
								if (vItem.size() != 5)
								{
									cout << string_format("'alpao pulse <act#> <amp> <PulseTime_ms>' Incorrect number of parameters (3 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								//sscanf(cmd_in, "%s %f", auxstr, &fChoice);
								//printf("Exposure time: %f\n", fChoice);
								int act_num = int(std::atoi(vItem[2].data()));
								float amp = float(std::atof(vItem[3].data()));
								float pulseTime_ms = float(std::atof(vItem[4].data()));

								cout << string_format("alpao pulse #: %d amp: %.6f PulseTime_ms: %.3e command issued.\n",
									act_num, amp, pulseTime_ms);
							}
							else if (compare_no_case(cmdStr,"poke") == 0)
							{
								if (vItem.size() != 4)
								{
									cout << string_format("'alpao poke <act#> <amp>' Incorrect number of parameters (2 required). Nothing done.");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								//sscanf(cmd_in, "%s %f", auxstr, &fChoice);
								//printf("Exposure time: %f\n", fChoice);
								int act_num = int(std::atoi(vItem[2].data()));
								float amp = float(std::atof(vItem[3].data()));

								cout << string_format("alpao poke #: %d amp: %.6f command issued.\n",
									act_num, amp);
							}
							else if (compare_no_case(cmdStr,"direct_dm_update") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("'alpao update_combined' incorrect number of parameters (1 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								if (compare_no_case(vItem[2], "on") == 0)
								{
									SetEvent(TP_DMcmd.hStopServer);
								}
								else if (compare_no_case(vItem[2], "off") == 0)
								{
									SetEvent(TP_DMcmd.hStartServer);
								}

								cout << string_format("'alpao direct_dm_update %s' command issued.\n", vItem[2].data());
							}
							else if (compare_no_case(cmdStr, "init_dm") == 0 && !bSimu_DM)
							{
								if (bSimu_DM == false)
								{

									if (bSimu_iXon == false)
									{
										if (dm.Check() != true)
										{
											system("pause");
											exit(-1);
										}
									}

									DM_nAct = dm.Get("NBOfActuator");
									cout << string_format("ALPAO DM with %d actuators initialized.\n\n",DM_nAct);

									if (DM_nAct != nAct)
									{
										cout << string_format("ALPAO DM detected # of actuators (%d) does not match the user input (%d). Quitting.\n",DM_nAct, nAct);
										system("pause");
										exit(-1);
									}

								}
								else
								{
									DM_nAct = nAct;
									cout << string_format("(SIMU) ALPAO DM with %d actuators initialized.\n\n",DM_nAct);
								}

								vDMSquare2DMmap.resize(DM_nAct, 0);
								vector<int> vSize;
								if (FITS.Load("R:\\ActuatorsMappingIndices.fits",vDMSquare2DMmap,vSize) != SUCCESS)
								{
									cout << string_format("'R:\\ActuatorsMappingIndices.fits' could not be loaded (check file existence). Quitting.\n");
									system("pause");
									exit(-1);
								}
		
								if (vSize[0] != DM_nAct)
								{
									cout << string_format("'R:\\ActuatorsMappingIndices.fits' mapping file size (%d) not matching the actual DM actuator number (%d). Quitting.\n",vSize[0],DM_nAct);
									system("pause");
									exit(-1);
								}

								// POSITION VALUE FOR EACH ACTUATOR (-1.0 <...< +1.0)
								vDM_cur_cmd.resize(DM_nAct, 0.f);
								vDM_static_cmd.resize(DM_nAct, 0.f);

								// DM COMMAND BUFFERS = 8x97 ARRAY
								vDM_cmd_buffers.resize(N_DM_COMMAND_MAPS);
								for (size_t k=0;k<N_DM_COMMAND_MAPS;++k)
								{
									vDM_cmd_buffers[k].resize(DM_nAct, 0.f);
								}

								// START THE DM UPDATE LOOP THROUGH COMMAND BUFFERS
								SetEvent(TP_DMcmd.hStartServer);

							}
							else if (compare_no_case(cmdStr,"load_flat") == 0)
							{
								if (vItem.size() < 3)
								{
									cout << string_format("ERROR: 'alpao load_flat \"<strFile>\"' Incorrect number of parameters (1 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								string strFile;
								concat_quoted_string_items(strFile, vItem, 2);

								if (MultiCODE::CheckFileExists(strFile) == false)
								{
									cout << string_format("ERROR: 'alpao load_flat \"<strFile>\"' could not be loaded (file does not exist). Nothing done.\n");
									break;
								}

								FITS.OpenFITS(strFile);
								string strVal;
								string comment;
								FITS.GetKeywordValue("ISFLAT",strVal,comment);
								if (std::atoi(strVal.data()) == 0)
								{
									cout << string_format("ERROR: 'alpao load_flat \"<strFile>\"' is not a flat file ('ISFLAT'=0). Nothing done.\n");
									break;
								}
								FITS.CloseFITS();

								vector<int> vSize;
								vector<float> vDM_flat;
								FITS.Load(strFile,vDM_flat,vSize);
								if (vSize[0] != TP_calcSH.nAct)
								{
									cout << string_format("ERROR: 'alpao load_flat': flat map must contain %d float values. Nothing done.\n", TP_calcSH.nAct);
									break;
								}

								std::copy(vDM_flat.begin(),vDM_flat.end(),vDM_cmd_buffers[0].begin());

								// SET MAP TO SHARED MEMORY BUFFER #0
								std::copy(vDM_flat.begin(),vDM_flat.end(),vSHM_DM_BUFFERS[0]->array.F);
								vSHM_DM_BUFFERS[0]->md[0].cnt0++;

								// update_DM_static_buffer DOES *NOT* UPDATE THE DM COMMAND MAP (IT JUST UPDATES THE STATIC MAP)
								update_DM_static_buffer();
								collapse_DM_buffers();

								cout << string_format("Alpao DM flat file '%s' loaded into slot #0.\n", strFile.data());
							}
							else if (compare_no_case(cmdStr,"load_map") == 0)
							{
								if (vItem.size() < 4)
								{
									cout << string_format("ERROR: 'alpao load_map \"<strFile>\" <channel#>' Incorrect number of parameters (2 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								string strFile;
								concat_quoted_string_items(strFile, vItem, 2);

								if (MultiCODE::CheckFileExists(strFile) == false)
								{
									cout << string_format("ERROR: 'alpao load_flat \"<strFile>\" <channel#>' could not be loaded (file does not exist). Nothing done.\n");
									break;
								}

								int channel_index = std::atoi(vItem[3].data());
								if ((channel_index < 0) || (channel_index > (N_DM_COMMAND_MAPS - 1)))
								{
									cout << string_format("ERROR: 'alpao load_flat \"<strFile>\" <channel#>' invalid channel number. Nothing done.\n");
									break;
								}

								FITS.OpenFITS(strFile);
								string strVal;
								string comment;
								FITS.GetKeywordValue("ISDMMAP",strVal,comment);
								if (std::atoi(strVal.data()) == 0)
								{
									cout << string_format("ERROR: 'alpao load_flat \"<strFile>\"' is not a valid DM map file ('ISDMMAP'=0). Nothing done.\n");
									break;
								}
								FITS.CloseFITS();

								vector<int> vSize;
								vector<float> vDM_map;
								FITS.Load(strFile,vDM_map,vSize);
								if (vSize[0] != TP_calcSH.nAct)
								{
									cout << string_format("ERROR: 'alpao load_map': DM map must contain %d float values. Nothing done.\n", TP_calcSH.nAct);
									break;
								}

								std::copy(vDM_map.begin(),vDM_map.end(),vDM_cmd_buffers[channel_index].begin());

								// SET MAP TO SHARED MEMORY BUFFER #0
								//for (size_t i=0;i<TP_calcSH.nAct;++i)
								//{
								//	vSHM_DM_BUFFERS[0]->array.F[vDMSquare2DMmap[i]] = vDM_flat[i];
								//}
								std::copy(vDM_map.begin(),vDM_map.end(),vSHM_DM_BUFFERS[channel_index]->array.F);
								vSHM_DM_BUFFERS[channel_index]->md[0].cnt0++;

								// update_DM_static_buffer DOES *NOT* UPDATE THE DM COMMAND MAP (IT JUST UPDATES THE STATIC MAP)
								update_DM_static_buffer();
								collapse_DM_buffers();

								cout << string_format("Alpao DM map file '%s' loaded into slot #%d.\n", strFile.data(), channel_index);
							}
							else if (compare_no_case(cmdStr, "shutdown_dm") == 0)
							{
								shutdown_dm();
							}
							else
							{
								cout << string_format("'%s' unknown command.\n", cmdStr.data());
							}

							break;		// TO EXIT THE while(true) STATEMENT
						}
						else if (compare_no_case(module,"algo") == 0)
						{
							break;		// TO EXIT THE while(true) STATEMENT
						}
						else if (compare_no_case(module,"system") == 0)
						{
							if (compare_no_case(cmdStr, "log") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								if (vItem[2].compare("ON") == 0)
								{
									SetEvent(TP_DataLogger.hStartLog);
									if (!TP_calcSH.isCal) {
										std::string timeDir;
										SYSTEMTIME myST;
										GetSystemTime(&myST);
										timeDir = string_format("%04d.%02d.%02dT%02dh%02dm%02d.%03ds", myST.wYear, myST.wMonth, myST.wDay,
											myST.wHour, myST.wMinute, myST.wSecond, myST.wMilliseconds);
									
										std::string strDirName = string_format("%s_%s_%s", TP_DataSaver.strPrefix.data(), timeDir.data(), "_NONE");
										TP_DataSaver.strLogDir = "C:\\Users\\lucas\\Documents\\STAGE\\Misc\\log\\" + strDirName + "\\";
										std::string strLogDirCalib = TP_DataSaver.strLogDir + "calib\\";
										_mkdir(TP_DataSaver.strLogDir.c_str());
										_mkdir(strLogDirCalib.c_str());
									}
								}
								else if (vItem[2].compare("OFF") == 0)
								{
									SetEvent(TP_DataLogger.hStopLog);
								}
								else if (vItem[2].compare("INIT") == 0)
								{
									SetEvent(TP_DataLogger.hInitializeLogBufSize);
								}
							}
							if (compare_no_case(cmdStr, "shutdown") == 0)
							{
								shutdown_dm();
							}
							else if (compare_no_case(cmdStr, "set_cmd_matrix") == 0)
							{
								if (vItem.size() != 2)
								{
									cout << string_format("Incorrect number of parameters (0 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								vector<int> vSize;
								FITS.Load("R:\\RRsvd.fits", vCmdMatrix, vSize);

								if ((vSize[0] != (TP_calcSH.nx * TP_calcSH.ny * 2)) || (vSize[1] != TP_calcSH.nModes_ZER))
								{
									cout << string_format("Invalid command matrix size (%dx%d) instead of %dx%d. Nothing done.\n",TP_calcSH.nx * TP_calcSH.ny * 2,TP_calcSH.nModes_ZER);
									break;		// TO EXIT THE while(true) STATEMENT
								}
							}
							//// -----------------------------------
							////       QUERY PARAMS
							//// -----------------------------------

							else if (compare_no_case(cmdStr, "query_params") == 0)
							{
								std::stringstream dest;
								dest << TP_calcSH.nAct;
								dest << ";";
								dest << TP_calcSH.nact_lin_x;
								dest << ";";
								dest << TP_calcSH.nact_lin_y;
								dest << ";";
								!bSimu_DM ? dest << "0" : dest << "1";
								dest << ";";
								!bSimu_iXon ? dest << "0" : dest << "1";
								std::string result = dest.str();
								std::sprintf(msgout, "%s", result.c_str());
								if (WriteFile(hPipeIn,
									msgout,
									std::strlen(msgout),   // = length of string + terminating '\0' !!!
									&dwWritten,
									NULL) == FALSE)
								{
									cout << string_format("ERROR: WriteFile(hPipeIn) failed!\n");
								}
								FlushFileBuffers(hPipeIn);
								
								break;
							}
							else if (compare_no_case(cmdStr, "simulation") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.\n");
									break;
								}
								if (vItem[2].compare("ON") == 0)
								{
									cout << "simulaton on" << endl;
									TP_iXon.bSimulation = true;
								}
								if (vItem[2].compare("OFF") == 0)
								{
									cout << "simulaton off" << endl;
									TP_iXon.bSimulation = false;
								}
							}
							else if (compare_no_case(cmdStr, "change_directory") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.\n");
									break;
								}
								else {
									SetEvent(TP_Simu_iXon.hStopLiveStream);
									SetEvent(TP_iXon.hStopAcquisition);
									Sleep(50);
									SimulationLoader::getInstance().emptyValues();
									SimulationLoader::getInstance().loadFullRepo(vItem[2] + '/');
									TP_Simu_iXon.maxImages = SimulationLoader::getInstance().getAllValues().size();
									TP_Simu_iXon.currentIndex = 0;
									SetEvent(TP_iXon.hStartAcquisition);
									SetEvent(TP_Simu_iXon.hStartLiveStream);
								}
							}
							else if (compare_no_case(cmdStr, "change_prefix") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.\n");
									break;
								}
								else {
									TP_DataSaver.strPrefix = vItem[2];
								}
							}	
							else if (compare_no_case(cmdStr, "set_ZER_gain") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								float gain = float(std::atof(vItem[2].data()));

								if ((gain < 0.f) || (gain > 10.f))
								{
									cout << string_format("Invalid gain value (0.0 < .. < 10.0). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								TP_calcSH.ZER_gain = gain;
								cout << string_format("'ZER set gain' command has been issued with value: %.5f.\n",gain);
							}
							else if (compare_no_case(cmdStr, "close_loop_ZER") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("'alpao update_combined' incorrect number of parameters (1 required). Nothing done.\n");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								if (compare_no_case(vItem[2], "on") == 0)
								{
									TP_calcSH.bEnableCloseLoop = true;
									cout << "'ZER closed loop ON' command has been issued.\n";
								}
								else if (compare_no_case(vItem[2], "off") == 0)
								{
									TP_calcSH.bEnableCloseLoop = false;
									//std::fill(vTT_drift.begin(),vTT_drift.end(),0.);
									//std::fill(vHO_drift.begin(),vHO_drift.end(),0.);
									cout << "'ZER closed loop OFF' command has been issued.\n";
								}

								cout << string_format("'alpao direct_dm_update %s' command issued.\n", vItem[2].data());

							}
							else if (compare_no_case(cmdStr,"calib_zer") == 0)
							{
								cout << string_format("System command '%s' issued.\n", cmdStr.data());

								SYSTEMTIME myST;
								GetSystemTime(&myST);
								std::string timeDir = string_format("%04d.%02d.%02dT%02dh%02dm%02d.%03ds", myST.wYear, myST.wMonth, myST.wDay,
									myST.wHour, myST.wMinute, myST.wSecond, myST.wMilliseconds);


								std::string strDirName = string_format("%s_%s_%s", TP_DataSaver.strPrefix.data(), timeDir.data(), "_ZER");
								TP_DataSaver.strLogDir = "C:\\Users\\lucas\\Documents\\STAGE\\Misc\\log\\" + strDirName + "\\";
								std::string strLogDirCalib = TP_DataSaver.strLogDir + "calib\\";
								_mkdir(TP_DataSaver.strLogDir.c_str());
								_mkdir(strLogDirCalib.c_str());


								if (TP_iXon.bAcqON)
								{
									//cconf->babort = true;
									SetEvent(TP_iXon.hStopAcquisition);
									if (WaitForSingleObject(TP_iXon.hAcquisitionDone,5000) == WAIT_TIMEOUT)
									{
										cout << string_format("TIMEOUT while waiting for acquisition routine to terminate\n");
									}

									//::AbortAcquisition(); // camera command
									TP_iXon.bAcqON = false;
									TP_iXon.bStreamON = false;

									cout << "iXon camera was streaming and has been stopped.\n";
								}

								if (TP_calcSH.bCalDataAvailable_ZER == true)
								{
									TP_calcSH.projMode = TP_calcSH.AOC_PROJ_MODE_ZER;
									TP_calcSH.bCalMode = true;

									// INITIALIZE DISPLAY BUFFERS
									TP_calcSH.display_image_count = TP_iXon.display_image_count;
									//TP_calcSH.vFlatSHData.resize(size_t(TP_calcSH.display_image_count*TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets));
									TP_calcSH.vModesDataHist.resize(TP_calcSH.display_image_count*TP_calcSH.MaxModes);

									vector<float> vtmp(TP_calcSH.display_image_count*TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets, 0.f);
									TP_calcSH.vCircBuf_FlatSHData.resize(TP_calcSH.nCircBuffers, vtmp);

									vCmdMatrix.resize(TP_calcSH.nModes_ZER * TP_calcSH.nx * TP_calcSH.ny * 2);

									for (size_t u=0;u<TP_calcSH.nModes_ZER;++u)
									{
										// SET DM POSITION BEFORE RECORDING SLOPES
										std::copy(TP_calcSH.vDM_Cmd_ZER.begin() + u*TP_calcSH.nAct, TP_calcSH.vDM_Cmd_ZER.begin() + (u+1)*TP_calcSH.nAct,
											vSHM_DM_BUFFERS[6]->array.F);
										vSHM_DM_BUFFERS[6]->md[0].cnt0++;

										if (true)
										{
											TP_calcSH.local_counter = 0;

											TP_calcSH.curCalMode = int(u);

											TP_iXon.bCalibration = true;
											TP_iXon.nCalFrames = TP_calcSH.nCalFrameCount;

											TP_iXon.bAcqON = true;
											TP_iXon.bStreamON = true;
											SetEvent(TP_iXon.hStartAcquisition);

											int waitResult = WaitForSingleObject(TP_iXon.hAcquisitionDone, 10000);

											switch (waitResult)
											{
											case WAIT_OBJECT_0:
												cout << string_format("Mode %d/%d acquired.\n", u + 1, TP_calcSH.nModes_ZER);
												SetEvent(TP_calcSH.hOut_SingleModeDone);
												break;
											case WAIT_TIMEOUT:
												cout << string_format("TIMEOUT.\n");
												break;
											}
										}

										Sleep(100);
									}

									std::fill(vSHM_DM_BUFFERS[6]->array.F,vSHM_DM_BUFFERS[6]->array.F + nAct, 0.f);
									vSHM_DM_BUFFERS[6]->md[0].cnt0++;

									int count = 0;
									for (size_t i = 0; i < TP_calcSH.nModes_ZER; ++i)
									{
										TP_calcSH.vDM_Modes_ZER_iNorm[i] = 0.f;
										for (size_t u = 0; u < size_t(TP_calcSH.nx * TP_calcSH.ny) * 2; ++u)
										{
											TP_calcSH.vDM_Modes_ZER_Slope[count] /= double(TP_calcSH.nCalFrameCount);
											TP_calcSH.vDM_Modes_ZER_iNorm[i] += TP_calcSH.vDM_Modes_ZER_Slope[count]*TP_calcSH.vDM_Modes_ZER_Slope[count];
											++count;
										}
										TP_calcSH.vDM_Modes_ZER_iNorm[i] = 1./TP_calcSH.vDM_Modes_ZER_iNorm[i];
									}

									TP_calcSH.isCal = true;

									FITS.Save(TP_calcSH.vDM_Modes_ZER_Slope,TP_calcSH.nx*TP_calcSH.nx*2,TP_calcSH.nModes_ZER,"R:\\","ZER_RecModes");
									FITS.Save(TP_calcSH.vDM_Modes_ZER_iNorm,TP_calcSH.nModes_ZER,1,"R:\\","ZER_iNormalizations");
									FITS.Save(TP_calcSH.vDM_Modes_ZER_Slope, TP_calcSH.nx * TP_calcSH.nx * 2, TP_calcSH.nModes_ZER, strLogDirCalib, "ZER_RecModes");
									FITS.Save(TP_calcSH.vDM_Modes_ZER_iNorm, TP_calcSH.nModes_ZER, 1, strLogDirCalib, "ZER_iNormalizations");

									TP_calcSH.bCalMode = false;
									TP_iXon.bCalibration = false;

									SetEvent(TP_EvtWatch.hOut_CalibDone);

									cout << string_format("ZERnike calibration done.\n");
								}
								else
								{
									cout << string_format("WARNING: ZERnike calibration data not available!\n");
								}

								break;
							}
							else if (compare_no_case(cmdStr, "query_log_dir") == 0)
							{
								std::sprintf(msgout, "%s", TP_DataSaver.strLogDir.c_str());
								if (WriteFile(hPipeIn,
									msgout,
									std::strlen(msgout),   // = length of string + terminating '\0' !!!
									&dwWritten,
									NULL) == FALSE)
								{
									cout << string_format("ERROR: WriteFile(hPipeIn) failed!\n");
								}
								FlushFileBuffers(hPipeIn);
							}
							else if (compare_no_case(cmdStr, "set_leakage_gain") == 0)
							{
								if (vItem.size() != 3)
								{
									cout << string_format("Incorrect number of parameters (1 required). Nothing done.");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								float _leakage_gain = float(std::atof(vItem[2].data()));

								if ((_leakage_gain < 0.f) || (_leakage_gain > 1.f))
								{
									cout << string_format("Invalid 'leakage_gain' value (0 <=...<= 1). Nothing done.");
									break;		// TO EXIT THE while(true) STATEMENT
								}

								leakage_gain = _leakage_gain;
								cout << string_format("'Set leakage gain' command has been issued with value: %.5f.\n",leakage_gain);
							}

							if (compare_no_case(cmdStr,"calib_zon") == 0)
							{
								if (TP_calcSH.bCalDataAvailable_ZON == true)
								{
									SYSTEMTIME myST;
									GetSystemTime(&myST);
									std::string timeDir = string_format("%04d.%02d.%02dT%02dh%02dm%02d.%03ds", myST.wYear, myST.wMonth, myST.wDay,
										myST.wHour, myST.wMinute, myST.wSecond, myST.wMilliseconds);


									std::string strDirName = string_format("%s_%s_%s", TP_DataSaver.strPrefix.data(), timeDir.data(), "_ZON");
									TP_DataSaver.strLogDir = "C:\\Users\\lucas\\Documents\\STAGE\\Misc\\log\\" + strDirName + "\\";
									std::string strLogDirCalib = TP_DataSaver.strLogDir + "calib\\";
									_mkdir(TP_DataSaver.strLogDir.c_str());
									_mkdir(strLogDirCalib.c_str());


									TP_calcSH.projMode = TP_calcSH.AOC_PROJ_MODE_ZON;
									TP_calcSH.bCalMode = true;

									cout << string_format("WARNING: ZON calibration not implemented (yet)!\n");
								}
								else
								{
									cout << string_format("WARNING: ZONnike calibration data not available!\n");
								}

								break;
							}

							break;
						}
						
						break;		// TO EXIT THE while(true) STATEMENT
					}
				}
				else
				{
					DisconnectNamedPipe(hPipeIn);
					CloseHandle(hPipeIn);
					hPipeIn = INVALID_HANDLE_VALUE;

					hPipeIn = CreateNamedPipe(TEXT(myfifoin),
					PIPE_ACCESS_DUPLEX,
					PIPE_TYPE_MESSAGE | PIPE_READMODE_BYTE | PIPE_WAIT,   // FILE_FLAG_FIRST_PIPE_INSTANCE is not needed but forces CreateNamedPipe(..) to fail if the pipe already exists...
					1,
					1024 * 64,
					1024 * 64,
					PipeWait,
					NULL);

					//Sleep(1000);
					break;
				}
			}
		}

		DisconnectNamedPipe(hPipeIn);
		DisconnectNamedPipe(hPipeOut);
	}


	printf("Shutting things down!\n");

	fflush(stdout);

	//free(imsize);
	//free(imarray);

	AOC_ShutdownCamera(&TP_iXon);

	//free(cconf);
	//cconf = NULL;				// free cam_config structure

	return;

}

int main(int argc, char *argv[])
{
	MCODE_Export_FITS FITS;
	//vector<string> vItem = {"module","cmd","\"C:\\Home\\Downloads Test\\Sub Folder1\\Test.fits\""};
	//string outstr;
	//concat_quoted_string_items(outstr,vItem,2);


	if (argc != 9)
	{
		return -1;
	}

	int _nThreads = std::atoi(argv[1]);
	int _nact_lin_x = std::atoi(argv[2]);
	int _nact_lin_y = std::atoi(argv[3]);
	int _nAct = std::atoi(argv[4]);
	int _DX = std::atoi(argv[5]);
	int _DY = std::atoi(argv[6]);
	bool _bSimu_iXon = std::atoi(argv[7]) == 0 ? false : true;
	bool _bSimu_DM = std::atoi(argv[8]) == 0 ? false : true;

	string tmpStr = string_format("\nAOC++ rev.%s\n=============\n\n",REVISION.data());
	cout << tmpStr;

	//tmpStr = string_format("Arguments: %d %d %d %d %d %s\n\n",_nThreads,_nact_lin_x,_nact_lin_y,_DX,_DY,_bSimulation == true ? "True" : "False");
	cout << string_format(
		"Arguments:\n"
		"  Number of threads (nThreads):                       %d\n"
		"  Linear number of DM 'x' actuators (nact_lin_x):     %d\n"
		"  Linear number of DM 'y' actuators (nact_lin_y):     %d\n"
		"  Number of DM actuators (nAct):                      %d\n"
		"  SH sensor 'x' number of pixels (DX):                %d\n"
		"  SH sensor 'y' number of pixels (DY):                %d\n"
		"  Simulation mode for iXon (bSimu_iXon):              %s\n"
		"  Simulation mode for DM (bSimu_DM):                  %s\n\n",
		_nThreads,_nact_lin_x,_nact_lin_y,_nAct,_DX,_DY,_bSimu_iXon == true ? "True" : "False",_bSimu_DM == true ? "True" : "False"
		);

	Log = new MCODE_LogManager("C:\\Users\\lucas\\Documents\\STAGE\\AOC_ServerLog.log",true);


	//string evtName = "_AOC_GENERAL_TRIGGER";
	//HANDLE testHandle = CreateEvent(NULL,FALSE,FALSE,evtName.data());

	bSimu_iXon = _bSimu_iXon;
	bSimu_DM = _bSimu_DM;
	nAct = _nAct;

	// START SHACK-HARTMANN CALCULUS MAIN THREAD
	TP_calcSH.nThreads = _nThreads;
	TP_calcSH.nact_lin_x = _nact_lin_x;
	TP_calcSH.nact_lin_y = _nact_lin_y;
	TP_calcSH.nAct = _nAct;
	TP_calcSH.nx = 10;
	TP_calcSH.ny = 10;
	TP_calcSH.nDM_buffers = N_DM_COMMAND_MAPS;
	TP_calcSH.DX = _DX;
	TP_calcSH.DY = _DY;
	TP_calcSH.sizex = _DX;
	TP_calcSH.sizey = _DY;
	TP_calcSH.SH_x0 = 0;
	TP_calcSH.SH_y0 = 0;
	TP_calcSH.fSH_dx = 12.8;
	TP_calcSH.fSH_dy = 12.8;
	TP_calcSH.MaxModes = 100;
	TP_calcSH.ZER_gain = 0.001;

	vDM_cmd.resize(TP_calcSH.nAct, 0.);
	vtmpLO.resize(TP_calcSH.nAct);
	vtmpHO.resize(TP_calcSH.nAct);

	vTT_drift.resize(2);
	vHO_drift.resize(TP_calcSH.MaxModes - 2);

	leakage_gain = 0.99f;

	// SET THE SH CELLS X BOUNDARIES
	for (int k=0;k<TP_calcSH.nact_lin_x+1;++k)
	{
		TP_calcSH.vSH_x.push_back(TP_calcSH.SH_x0 + int(k * TP_calcSH.fSH_dx + 0.5));
	}

	// SET THE SH CELLS Y BOUNDARIES
	for (int k=0;k<TP_calcSH.nact_lin_y+1;++k)
	{
		TP_calcSH.vSH_y.push_back(TP_calcSH.SH_y0 + int(k * TP_calcSH.fSH_dy + 0.5));
	}
	
	TP_calcSH.vSH_ValidCell.resize(TP_calcSH.nx*TP_calcSH.ny, 1);
	TP_calcSH.vSH_ValidCell[0] = 0;
	TP_calcSH.vSH_ValidCell[1] = 0;
	TP_calcSH.vSH_ValidCell[2] = 0;
	TP_calcSH.vSH_ValidCell[7] = 0;
	TP_calcSH.vSH_ValidCell[8] = 0;
	TP_calcSH.vSH_ValidCell[9] = 0;
	TP_calcSH.vSH_ValidCell[10] = 0;
	TP_calcSH.vSH_ValidCell[11] = 0;
	TP_calcSH.vSH_ValidCell[18] = 0;
	TP_calcSH.vSH_ValidCell[19] = 0;
	TP_calcSH.vSH_ValidCell[20] = 0;
	TP_calcSH.vSH_ValidCell[29] = 0;
	TP_calcSH.vSH_ValidCell[70] = 0;
	TP_calcSH.vSH_ValidCell[79] = 0;
	TP_calcSH.vSH_ValidCell[80] = 0;
	TP_calcSH.vSH_ValidCell[81] = 0;
	TP_calcSH.vSH_ValidCell[88] = 0;
	TP_calcSH.vSH_ValidCell[89] = 0;
	TP_calcSH.vSH_ValidCell[90] = 0;
	TP_calcSH.vSH_ValidCell[91] = 0;
	TP_calcSH.vSH_ValidCell[92] = 0;
	TP_calcSH.vSH_ValidCell[97] = 0;
	TP_calcSH.vSH_ValidCell[98] = 0;
	TP_calcSH.vSH_ValidCell[99] = 0;
	// 76 VALID SUBAPERTURES

	if (StartThread_calcSH(_nThreads, TP_calcSH.DX, TP_calcSH.DY, TP_calcSH.nx, TP_calcSH.ny) != SUCCESS)
	{
		cout << string_format("Failed to start thread 'calcSH'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	if (StartThread_DMcmd() != SUCCESS)
	{
		cout << string_format("Failed to start thread 'DMcmd'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	int nCells = TP_calcSH.nx * TP_calcSH.ny;
	int nsub = nCells / _nThreads;
	//vector<mcUINT32> vInd_i((_nThreads-1) * nsub);
	//vector<mcUINT32> vInd_j((_nThreads-1) * nsub);
	//vector<vector<int>> vInd_i(_nThreads - 1);
	//vector<vector<int>> vInd_j(_nThreads - 1);
	TP_calcSH.vInd_i.resize(size_t(_nThreads - 1));
	TP_calcSH.vInd_j.resize(size_t(_nThreads - 1));
	mcUINT32 cellCounter = 0;

	TP_calcSH.vSH_cent_x.resize(nCells,0.f);
	TP_calcSH.vSH_cent_y.resize(nCells,0.f);
	TP_calcSH.vSH_slope_x.resize(nCells,0.f);
	TP_calcSH.vSH_slope_y.resize(nCells,0.f);
	TP_calcSH.vSH_phot.resize(nCells,0.f);
	TP_calcSH.vThresh.resize(nCells,0.f);
	TP_calcSH.vEnabled.resize(nCells,1);

	// SET THE DEFAULT X REFERENCE POSITIONS (CENTER OF CELL)
	for (int j = 0; j < TP_calcSH.ny; ++j)
	{
		for (int i = 0; i < TP_calcSH.nx; ++i)
		{
			TP_calcSH.vSH_ref_x.push_back(0.5 * (TP_calcSH.vSH_x[i+1] - TP_calcSH.vSH_x[i]));
			TP_calcSH.vSH_ref_y.push_back(0.5 * (TP_calcSH.vSH_y[j+1] - TP_calcSH.vSH_y[j]));
		}
	}

	FITS.Save(TP_calcSH.vSH_ref_x,int(TP_calcSH.vSH_ref_x.size()),1,"R:\\","default_ref_pos_x");
	FITS.Save(TP_calcSH.vSH_ref_y,int(TP_calcSH.vSH_ref_y.size()),1,"R:\\","default_ref_pos_y");

	TP_calcSH.logVerbose = 0;

	for (int k=0;k<_nThreads - 1;++k)
	{
		TP_calcSH.vInd_i[k].resize(nsub);
		TP_calcSH.vInd_j[k].resize(nsub);
		for (size_t u=0;u<nsub;++u)
		{
			TP_calcSH.vInd_i[k][u] = cellCounter % TP_calcSH.nx;
			TP_calcSH.vInd_j[k][u] = cellCounter / TP_calcSH.nx;
			TP_calcSH.vThresh[cellCounter] = 0.f;
			TP_calcSH.vEnabled[cellCounter] = 1;
			++cellCounter;
		}
		//set_thread_params(k,nsub,&TP_calcSH.vSH_x[0],&TP_calcSH.vSH_y[0],&vInd_i[k][0],&vInd_j[k][0],&TP_calcSH.vSH_cent_x[0],&TP_calcSH.vSH_cent_y[0],&TP_calcSH.vSH_phot[0]);
		set_thread_params(k,nsub,TP_calcSH.vSH_x,TP_calcSH.vSH_y,TP_calcSH.vInd_i[k],TP_calcSH.vInd_j[k],
			TP_calcSH.vSH_cent_x,TP_calcSH.vSH_cent_y,TP_calcSH.vSH_slope_x,TP_calcSH.vSH_slope_y,TP_calcSH.vSH_phot,
			TP_calcSH.vSH_ref_x,TP_calcSH.vSH_ref_y,TP_calcSH.vThresh,TP_calcSH.vEnabled);
	}

	int rem_sub = nCells - nsub * (TP_calcSH.nThreads - 1);
	TP_calcSH.vInd_i2.resize(rem_sub);
	TP_calcSH.vInd_j2.resize(rem_sub);
	//TP_calcSH.vThresh2.resize(rem_sub);
	for (size_t u=0;u<rem_sub;++u)
	{
		TP_calcSH.vInd_i2[u] = cellCounter % TP_calcSH.nx;
		TP_calcSH.vInd_j2[u] = cellCounter / TP_calcSH.nx;
		//TP_calcSH.vThresh2[u] = 0.f;
		TP_calcSH.vThresh[cellCounter] = 0.f;
		++cellCounter;
	}
	//set_thread_params(TP_calcSH.nThreads - 1,rem_sub,&TP_calcSH.vSH_x[0],&TP_calcSH.vSH_y[0],&vInd_i2[0],&vInd_j2[0],&TP_calcSH.vSH_cent_x[0],&TP_calcSH.vSH_cent_y[0],&TP_calcSH.vSH_phot[0]);
	set_thread_params(TP_calcSH.nThreads - 1,rem_sub,TP_calcSH.vSH_x,TP_calcSH.vSH_y,TP_calcSH.vInd_i2,TP_calcSH.vInd_j2,
		TP_calcSH.vSH_cent_x,TP_calcSH.vSH_cent_y,TP_calcSH.vSH_slope_x,TP_calcSH.vSH_slope_y,TP_calcSH.vSH_phot,
		TP_calcSH.vSH_ref_x,TP_calcSH.vSH_ref_y,TP_calcSH.vThresh,TP_calcSH.vEnabled);


	//TP_EvtWatch.mapHandle.insert(std::pair<string, HANDLE>(evtName, testHandle));

	// START DATA LOGGER THREAD
	if (StartThread_DataLogger() != SUCCESS)
	{
		cout << string_format("Failed to start thread 'DataLogger'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	// START DATA SAVER THREAD
	if (StartThread_DataSaver() != SUCCESS)
	{
		cout << string_format("Failed to start thread 'DataSaver'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	// START ANDOR IXON CAMERA MAIN THREAD
	if (StartThread_EvtWatch() != SUCCESS)
	{
		cout << string_format("Failed to start thread 'EvtWatch'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	// START ANDOR IXON CAMERA MAIN THREAD
	if (StartThread_iXon() != SUCCESS)
	{
		cout << string_format("Failed to start thread 'iXon'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	// START SIMULATION IXON MAIN THREAD
	if (StartThread_Simu_iXon() != SUCCESS)
	{
		cout << string_format("Failed to start thread 'Simu iXon'! Quitting...") << endl;
		system("pause");
		exit(-1);
	}


	TP_iXon.bSimulation = _bSimu_iXon;

	// CamID = -1 IS ALSO THE DEFAULT. IF ONLY ONE CAMERA IS CONNECTED, CamID = 0 (SHOULD BE)
	if (_bSimu_iXon == true)
	{
		TP_iXon.camID = -1;
	}
	else
	{
		TP_iXon.camID = 0;
	}

	if (AOC_InitializeCamera(&TP_iXon) != SUCCESS)
	{
		cout << string_format("Failed to initialize iXon camera! Quitting...") << endl;
		system("pause");
		exit(-1);
	}

	// INITIALIZE DEFAULT HW PARAMETERS
	if (AOC_SetHWConf(&TP_iXon) != SUCCESS)
	{
		cout << string_format("Failed to initialize hardware configuration! Quitting...") << endl;
		system("pause");
		exit(-1);
	}


	// ------------------------------------------
	//          setup shared memory 
	// ------------------------------------------
	long naxis = 2;          // number of axes
	uint8_t atype;           // data type code
	uint32_t *imsize;        // image size
	int shared;              // set to 1 if image is in shared memory
	int NBkw;                // number of keywords supported

	vCollapseStaticInd.push_back(0);		// STATIC DM FLAT
	vCollapseStaticInd.push_back(5);		// STATIC 'POKER' COMMANDS
	vCollapseStaticInd.push_back(4);		// STATIC 'ZALAMANO' COMMANDS

	bool bFound = false;

	for (size_t k=0;k<TP_calcSH.nDM_buffers;++k)
	{
		bFound = false;

		for (size_t u=0;u<vCollapseStaticInd.size();++u)
		{
			if (vCollapseStaticInd[u] == k)
			{
				bFound = true;
				break;
			}
		}

		if (bFound == true)
		{
			continue;
		}

		vCollapseDynamicInd.push_back(k);
	}

	// AOC DM STATE: IT WILL CONTAIN THE COLLAPSED DATA FROM ALL CHANNELS
	naxis		= 2;                // 2D image
	SHM_DM_state	= (IMAGE*) malloc(sizeof(IMAGE));
	imsize		= (uint32_t *) malloc(naxis * sizeof(uint32_t));
	imsize[0]	= TP_calcSH.nAct;       // image dimensions
	imsize[1]	= 1;       // image dimensions
	atype		= _DATATYPE_FLOAT; // camera SDK writes "int"
	shared		= 1;               // image will be in shared memory
	NBkw		= 2;              // allocate space for 2 keywords
	ImageCreate(&SHM_DM_state[0], "aoc_dm", naxis, imsize, atype, shared, NBkw);

	// SH IXON IMAGE DISPLAY
	naxis		= 2;                // 2D image
	SHM_imarray	= (IMAGE*) malloc(sizeof(IMAGE));
	imsize		= (uint32_t *) malloc(naxis * sizeof(uint32_t));
	imsize[0]	= TP_calcSH.DX;       // image dimensions
	imsize[1]	= TP_calcSH.DY;       // image dimensions
	atype		= _DATATYPE_INT32; // camera SDK writes "int"
	shared		= 1;               // image will be in shared memory
	NBkw		= 2;              // allocate space for 2 keywords
	ImageCreate(&SHM_imarray[0], "ixon", naxis, imsize, atype, shared, NBkw);
	//SHM_imarray->kw[0].type = 'L';
	//strcpy(&SHM_imarray->kw[0].name[0],"exptime_microsecond");
	//SHM_imarray->kw[0].value.numl = int(TP_iXon.expTime_s * 1000000 + 0.5f);

	// ALL SH nx*ny*n SUBAPERTURE SLOPES DURING 'n' FRAMES
	TP_calcSH.nDataSets = 5;
	naxis		= 2;                // 2D image
	SHM_Slopes	= (IMAGE*) malloc(sizeof(IMAGE));
	imsize[0]	= TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets;		// image dimensions. x2 FOR X AND Y CENTROIDS, x2 FOR SLOPES, ILLUMINATION
	imsize[1]	= MAX_BUFFER_HISTORY_LENGTH;								// MAXIMUM HISTORY BUFFER SIZE
	atype		= _DATATYPE_FLOAT;	// camera SDK writes "int"
	shared		= 1;				// image will be in shared memory
	NBkw		= 2;				// allocate space for 2 keywords
	ImageCreate(&SHM_Slopes[0], "slope_history", naxis, imsize, atype, shared, NBkw);

	// ALLOCATE MEMORY FOR TIMESTAMP BUFFER
	vTimeStamp.resize(MAX_BUFFER_HISTORY_LENGTH, long long(0));

	vector<float> vtmp(TP_iXon.display_image_count*TP_calcSH.nx*TP_calcSH.ny*TP_calcSH.nDataSets, 0.f);
	vector<double> vtimestamptmp(TP_iXon.display_image_count * 50, 0);
	vector<float> vdmtmp(TP_iXon.display_image_count * 97, 0);
	TP_calcSH.nCircBuffers = 2;
	TP_calcSH.vCircBuf_FlatSHData.resize(TP_calcSH.nCircBuffers, vtmp);
	TP_calcSH.vCircBuf_Timestamps.resize(TP_calcSH.nCircBuffers, vtimestamptmp);
	TP_calcSH.vCircBuf_DM_cmd.resize(TP_calcSH.nCircBuffers, vdmtmp);
	TP_calcSH.circBufCounter = 0;

	// 
	naxis		= 2;                // 2D image
	SHM_RefCent	= (IMAGE*) malloc(sizeof(IMAGE));
	imsize[0]	= TP_calcSH.nx*TP_calcSH.ny*2;       // image dimensions. x2 FOR X AND Y SLOPES
	imsize[1]	= 1;				// image dimensions
	atype		= _DATATYPE_FLOAT;	// camera SDK writes "int"
	shared		= 1;				// image will be in shared memory
	NBkw		= 2;				// allocate space for 2 keywords
	ImageCreate(&SHM_RefCent[0], "reference_centroids", naxis, imsize, atype, shared, NBkw);

	// SET SHM REFERENCE POSITIONS TO PREVIOUSLY INITIALIZED (DEFAULT) VALUES
	std::copy(TP_calcSH.vSH_ref_x.begin(),TP_calcSH.vSH_ref_x.end(), (float*)SHM_RefCent[0].array.SI32);
	std::copy(TP_calcSH.vSH_ref_y.begin(),TP_calcSH.vSH_ref_y.end(), (float*)SHM_RefCent[0].array.SI32 + TP_calcSH.nx*TP_calcSH.ny);


	// 
	naxis		= 2;                // 2D image
	SHM_DM_slope_mask	= (IMAGE*) malloc(sizeof(IMAGE));
	imsize[0]	= TP_calcSH.nx*TP_calcSH.ny;       // image dimensions.
	imsize[1]	= 1;				// image dimensions
	atype		= _DATATYPE_INT32;	// camera SDK writes "int"
	shared		= 1;				// image will be in shared memory
	NBkw		= 2;				// allocate space for 2 keywords
	ImageCreate(&SHM_DM_slope_mask[0], "pupil_slope_mask", naxis, imsize, atype, shared, NBkw);


	// SET THE SHARED MEMORY BUFFER FOR COMBINED CHANNELS
	vSHM_DM_BUFFERS.resize(N_DM_COMMAND_MAPS);
	for (size_t k=0;k<N_DM_COMMAND_MAPS;++k)
	{
		naxis = 2;							// 2D image
		vSHM_DM_BUFFERS[k] = (IMAGE*) malloc(sizeof(IMAGE));
		imsize[0]	= TP_calcSH.nAct;       // image dimensions.
		imsize[1]	= 1;					// image dimensions
		atype		= _DATATYPE_FLOAT;		// camera SDK writes "int"
		shared		= 1;					// image will be in shared memory
		NBkw		= 30;					// allocate space for 30 keywords
		ImageCreate(&vSHM_DM_BUFFERS[k][0], string_format("DM_Buffer_%02d",k).data(), naxis, imsize, atype, shared, NBkw);
	}
	// TIP-TILT AMPLITUDE ACCUMULATED VALUES (BUFER #2)
	strcpy(&vSHM_DM_BUFFERS[2]->kw[0].name[0],"DM_Abs_TTx");
	strcpy(&vSHM_DM_BUFFERS[2]->kw[1].name[0],"DM_Abs_TTy");
	vSHM_DM_BUFFERS[2]->kw[0].type = 'D';
	vSHM_DM_BUFFERS[2]->kw[1].type = 'D';

	// HIGH ORDER POLYNOMIALS AMPLITUDE ACCUMULATED VALUES (BUFER #3)
	for (size_t k = 0; k < 30; ++k)
	{
		strcpy(&vSHM_DM_BUFFERS[3]->kw[k].name[0], string_format("DM_Abs_ZER_%02d",k).data());
		vSHM_DM_BUFFERS[3]->kw[k].type = 'D';
	}

	// SET THE SHARED MEMORY BUFFER FOR MODE-*RELATED INFORMATION
	naxis = 2;							// 2D image
	SHM_ModesInfo = (IMAGE*) malloc(sizeof(IMAGE));
	imsize[0]	= TP_calcSH.MaxModes;	// image dimensions.
	imsize[1]	= 500;					// MAXIMUM HISTORY BUFFER SIZE
	atype		= _DATATYPE_FLOAT;		// camera SDK writes "int"
	shared		= 1;					// image will be in shared memory
	NBkw		= 2;					// allocate space for 2 keywords
	ImageCreate(&SHM_ModesInfo[0], "modes_info", naxis, imsize, atype, shared, NBkw);


	// SET THE SHARED MEMORY BUFFER FOR TIMESTAMP
	naxis = 2;							// 2D image
	SHM_Timestamps = (IMAGE*)malloc(sizeof(IMAGE));
	imsize[0] = 50 * TP_calcSH.nDataSets;	// image dimensions.
	imsize[1] = 500;					// MAXIMUM HISTORY BUFFER SIZE
	atype = _DATATYPE_DOUBLE;		// camera SDK writes "int"
	shared = 1;					// image will be in shared memory
	NBkw = 2;					// allocate space for 2 keywords
	ImageCreate(&SHM_Timestamps[0], "timestamps", naxis, imsize, atype, shared, NBkw);



	//SetEvent(TP_DataLogger.hStartLog);

	// START THE MAIN LOOP (COMMAND INTERFACE)
	main_loop();


	SetEvent(TP_DataLogger.hStopLog);


	if (hThread_DMcmd != INVALID_HANDLE_VALUE)
	{
		TerminateThread_DMcmd();
	}

	if (hThread_iXon != INVALID_HANDLE_VALUE)
	{
		TerminateThread_iXon();
	}

	Log->Stop();
	delete Log;

	system("pause");

	return 0;
}