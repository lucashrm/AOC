// DEV_CIAO.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"

#include <MultiCODE.h>
#include <MCODE.h>
#include <MCODE_Export.h>
#include <MCODE_Export_FITS.h>
#include <MCODE_SciLib.h>
#include <MCODE_StdFunc.h>
#include <MCODE_SystemUtilities.h>

using namespace std;
using namespace MultiCODE;

#if defined(WIN64)
	#if defined(_DEBUG)
		#pragma comment(lib,"MCODE_Export_x64_D.lib")
		#pragma comment(lib,"MCODE_Export_FITS_x64_D.lib")
		#pragma comment(lib,"MCODE_Module_x64_D.lib")
		#pragma comment(lib,"MCODE_PThread_x64_D.lib")
		#pragma comment(lib,"MCODE_SciLib_x64_D.lib")
		#pragma comment(lib,"MCODE_StdFunc_x64_D.lib")
		#pragma comment(lib,"MCODE_SystemUtilities_x64_D.lib")
	#else
		#pragma comment(linker,"/NODEFAULTLIB:libcmt.lib")
		#pragma comment(lib,"MCODE_Export_x64.lib")
		#pragma comment(lib,"MCODE_Export_FITS_x64.lib")
		#pragma comment(lib,"MCODE_Module_x64.lib")
		#pragma comment(lib,"MCODE_PThread_x64.lib")
		#pragma comment(lib,"MCODE_SciLib_x64.lib")
		#pragma comment(lib,"MCODE_StdFunc_x64.lib")
		#pragma comment(lib,"MCODE_SystemUtilities_x64.lib")
	#endif
#elif defined(WIN32)
	#if defined(_DEBUG)
		#pragma comment(lib,"MCODE_Export_D.lib")
		#pragma comment(lib,"MCODE_Export_FITS_D.lib")
		#pragma comment(lib,"MCODE_Module_D.lib")
		#pragma comment(lib,"MCODE_PThread_D.lib")
		#pragma comment(lib,"MCODE_SciLib_D.lib")
		#pragma comment(lib,"MCODE_StdFunc_D.lib")
		#pragma comment(lib,"MCODE_SystemUtilities_D.lib")
	#else
		#pragma comment(linker,"/NODEFAULTLIB:libcmt.lib")
		#pragma comment(lib,"MCODE_Export.lib")
		#pragma comment(lib,"MCODE_Export_FITS.lib")
		#pragma comment(lib,"MCODE_Module.lib")
		#pragma comment(lib,"MCODE_PThread.lib")
		#pragma comment(lib,"MCODE_SciLib.lib")
		#pragma comment(lib,"MCODE_StdFunc.lib")
		#pragma comment(lib,"MCODE_SystemUtilities.lib")
	#endif
#endif

#if defined(WINDOWS)
	#pragma comment(lib,"pthreadVC2.lib")
#endif

#include <pthread.h>

//=================================
#define		SX		128
#define		SY		128

MCODE_SciLib sci;

typedef struct _ThreadParams_WFScalc
{
	int nx;						// TOTAL NUMBER OF SUBAPERTURES IN THE 'X' DIRECTION
	int ny;						// TOTAL NUMBER OF SUBAPERTURES IN THE 'Y' DIRECTION
	int sx;						// X SIZE OF THE SUB APERTURE [PIXELS]
	int sy;						// Y SIZE OF THE SUB APERTURE [PIXELS]
	int thresh;					// SUBAPERTURE INTENSITY THRESHOLD FOR CENTROID CALCULATION
	vector<float> vx;			// LINEAR VECTOR USED FOR Y CENTROID COMPUTATION
	vector<float> vy;			// LINEAR VECTOR USED FOR X CENTROID COMPUTATION
	vector<int> vSubArr;		// VECTOR WITH SUBAPERTURE DATA (MUST BE RESIZED PROPERLY BEFORE BEING USED)
	//sci.LinearVector(vx,SX);
	//sci.LinearVector(vy,SY);

} ThreadParams_WFScalc;

int calc_SH_data(vector<int>& v_xx, vector<int>& v_yy, vector<int>& vImg, vector<float>& SH_xtmp, vector<float>& SH_ytmp)
{
//def calc_SH_data(self, data=None, ref=True):
//    xx, yy   = self.SH_xx , self.SH_yy        
//    ncx, ncy = xx.size - 1, yy.size - 1
	int nx = int(v_xx.size());
	int ny = int(v_yy.size());
	int ncx = nx - 1;
	int ncy = ny - 1;

//    if data is None:
//        self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True)
//    else:
//        self.live_img = data

//    self.im_cnt = self.shm_im.get_counter()  # image counter
//    bckgd = self.live_img.mean()

	mcUINT32 bkg_ui = 0;
	for (size_t k=0;k<SX*SY;++k)
	{
		bkg_ui += vImg[k];
	}
	mcFLOAT32 mean_bkg = float(bkg_ui)/float(SX*SY);

//    #self.live_img[self.live_img <= self.vmin] = self.vmin

//    for j in xrange(ncy):
//        y0, y1 = int(np.round(yy[j])), int(np.round(yy[j+1]))

//        for i in xrange(ncx):
//            x0, x1 = int(np.round(xx[i])), int(np.round(xx[i+1]))

//            sub_arr           = self.live_img[y0:y1,x0:x1]
//            self.SH_phot[j,i] = sub_arr.max()# - bckgd# max()

//            sub_arr[sub_arr < self.threshold] = self.threshold
//            #if self.SH_phot[j,i] > self.threshold: # original line
//            if self.SH_phot[j,i] > 1.3 * bckgd: # 30 % above background?
//                (yc, xc) = centroid_position_0(sub_arr)
//            else:
//                (yc, xc) = (self.SH_xref[j,i], self.SH_yref[j,i])

//            self.SH_xtmp[j,i] = xc
//            self.SH_ytmp[j,i] = yc

	int x0;
	int y0;
	int x1;
	int y1;

	vector<int> sub_arr;

	double maxval;
	double thresh = double(1.3 * mean_bkg);

	//vector<float> vx;
	//vector<float> vy;
	//sci.LinearVector(vx,SX);
	//sci.LinearVector(vy,SY);

	SH_xtmp.resize(nx*ny);
	SH_ytmp.resize(nx*ny);

	int xs;
	int ys;

	double xc;
	double yc;

	LARGE_INTEGER t0;
	LARGE_INTEGER t1;

	LARGE_INTEGER liFreq;
	QueryPerformanceFrequency(&liFreq);
	QueryPerformanceCounter(&t0);

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

			sci.SubArray(sub_arr,SX,SY,vImg,x0,y0,xs,ys);

			maxval = sci.IDL_MAX(sub_arr);

			int crd = i + j*nx;
			if (maxval > thresh)
			{
				// SUBTRACT MIN
				for (int k=0;k<xs*ys;++k)
				{
					if (sub_arr[k] < thresh) sub_arr[k] = 0;
				}

				// COMPUTE CENTROID
				sci.Centroid(sub_arr,xs,ys,xc,yc);
			}
			else
			{
				// SET TO REF POSITION
				xc = 0.;
				yc = 0.;
			}

			SH_xtmp[crd] = float(x0 + xc);
			SH_ytmp[crd] = float(y0 + yc);
		}
	}

	QueryPerformanceCounter(&t1);

	cout << "Execution time: " << double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart) << endl;

//    if ref is True:
//        self.SH_xref = self.SH_xtmp.copy()
//        self.SH_yref = self.SH_ytmp.copy()

//    self.SH_xslp = self.SH_xtmp - self.SH_xref
//    self.SH_yslp = self.SH_ytmp - self.SH_yref

//    self.SH_phot -= self.threshold
//    self.SH_xslp[self.SH_phot < 0] = 0.0
//    self.SH_yslp[self.SH_phot < 0] = 0.0

//    self.SH_comb[0] = self.SH_xslp
//    self.SH_comb[1] = self.SH_yslp
//    self.SH_comb[2] = self.SH_phot / self.SH_phot.max()

//    self.shm_comb.set_data(self.SH_comb)
//    #self.shm_xslp.set_data(self.SH_xslp)
//    #self.shm_yslp.set_data(self.SH_yslp)
//    self.shm_phot_inst.set_data(self.SH_phot)
//    
//    # here is information about the tip-tilt in pixels!
//    # weighted mean version!
//    self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
//    self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)
//    #self.ttx_mean = np.sum(self.SH_xslp * self.SH_phot) / np.sum(self.SH_phot)
//    #self.tty_mean = np.sum(self.SH_yslp * self.SH_phot) / np.sum(self.SH_phot)

//    # original version below
//    #self.ttx_mean = np.median(self.SH_xslp[self.SH_phot > self.threshold])
//    #self.tty_mean = np.median(self.SH_yslp[self.SH_phot > self.threshold])

//    self.update_log()

	return SUCCESS;
}

int calc_SH_data(int sizex, int sizey, int *vImg, int nx, int ny, int *v_xx, int *v_yy, float *SH_xtmp, float *SH_ytmp)
{
	//LARGE_INTEGER t0;
	//LARGE_INTEGER t1;

	//LARGE_INTEGER liFreq;
	//QueryPerformanceFrequency(&liFreq);
	int ncx = nx - 1;
	int ncy = ny - 1;

	//QueryPerformanceCounter(&t0);

	int bkg_ui = 0;
	for (int k=0;k<sizex*sizey;++k)
	{
		bkg_ui += vImg[k];
	}
	float mean_bkg = float(bkg_ui)/float(sizex*sizey);

	int x0;
	int y0;
	int x1;
	int y1;

	vector<int> sub_arr;

	int maxval;
	double thresh = double(1.3 * mean_bkg);

	//vector<float> vx;
	//vector<float> vy;
	//sci.LinearVector(vx,SX);
	//sci.LinearVector(vy,SY);

	//SH_xtmp.resize(nx*ny);
	//SH_ytmp.resize(nx*ny);

	int xs;
	int ys;

	double xc;
	double yc;

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

			// EXTRACT THE SUBPUPIL DATA FROM ORIGINAL IMAGE
			//sci.SubArray(sub_arr,SX,SY,vImg,x0,y0,xs,ys);
			sub_arr.resize(xs*ys);

			for (size_t v=0;v<size_t(ys);++v)
			{
				for (size_t u=0;u<size_t(xs);++u)
				{
					sub_arr[u + v*xs] = vImg[x0+u + (y0+v)*sizex];
				}
			}
			
			// FIND THE MAXIMUM VALUE IN THE SUBPUPIL
			//maxval = sci.IDL_MAX(sub_arr);
			maxval = -INT_MAX;
			vector<int>::iterator ptr = sub_arr.begin();
			for (int k=0;k<xs*ys;++k, ++ptr)
			{
				if (*ptr > maxval) maxval = *ptr;
			}

			int crd = i + j*ncx;

			float xc = 0.f;
			float yc = 0.f;
			float sum = 0.f;

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
				//sci.Centroid(sub_arr,xs,ys,xc,yc);
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

			SH_xtmp[crd] = float(x0 + xc);
			SH_ytmp[crd] = float(y0 + yc);
		}
	}

	//QueryPerformanceCounter(&t1);

	//cout << "Execution time overall: " << double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart) << endl;

	return SUCCESS;
}


int calc_SH_data_mt(ThreadParams_WFScalc *param, const int i_xx, const int i_yy, const vector<int>& vImg, vector<float>& SH_xtmp, vector<float>& SH_ytmp)
{
//    if data is None:
//        self.live_img = self.shm_im.get_data(check=self.im_cnt,reform=True)
//    else:
//        self.live_img = data

//    self.im_cnt = self.shm_im.get_counter()  # image counter
//    bckgd = self.live_img.mean()

	mcUINT32 bkg_ui = 0;
	int offset = param->nx - param->sx;

	vector<int>::const_iterator data_ptr = vImg.begin() + (i_xx + i_yy*param->nx);
	vector<int>::iterator subArr_ptr = param->vSubArr.begin();
	
	for (int v=0;v<param->sy;++v)
	{
		for (int u=0;u<param->sx;++u)
		{
			*(subArr_ptr++) = *(data_ptr);
			bkg_ui += *(data_ptr++);
		}
		data_ptr += offset;
	}
	mcFLOAT32 mean_bkg = float(bkg_ui)/float(param->sx*param->sy);

//    #self.live_img[self.live_img <= self.vmin] = self.vmin

//    for j in xrange(ncy):
//        y0, y1 = int(np.round(yy[j])), int(np.round(yy[j+1]))

//        for i in xrange(ncx):
//            x0, x1 = int(np.round(xx[i])), int(np.round(xx[i+1]))

//            sub_arr           = self.live_img[y0:y1,x0:x1]
//            self.SH_phot[j,i] = sub_arr.max()# - bckgd# max()

//            sub_arr[sub_arr < self.threshold] = self.threshold
//            #if self.SH_phot[j,i] > self.threshold: # original line
//            if self.SH_phot[j,i] > 1.3 * bckgd: # 30 % above background?
//                (yc, xc) = centroid_position_0(sub_arr)
//            else:
//                (yc, xc) = (self.SH_xref[j,i], self.SH_yref[j,i])

//            self.SH_xtmp[j,i] = xc
//            self.SH_ytmp[j,i] = yc

	int crd = i_xx + i_yy*param->nx;

	double maxval;
	double thresh = double(1.3 * mean_bkg);

	//vector<float> vx;
	//vector<float> vy;
	//sci.LinearVector(vx,SX);
	//sci.LinearVector(vy,SY);

	SH_xtmp.resize(param->nx*param->ny);
	SH_ytmp.resize(param->nx*param->ny);

	int xs;
	int ys;

	double xc;
	double yc;

	//LARGE_INTEGER t0;
	//LARGE_INTEGER t1;

	//LARGE_INTEGER liFreq;
	//QueryPerformanceFrequency(&liFreq);
	//QueryPerformanceCounter(&t0);

	{
		vector<int>::iterator ptr = param->vSubArr.begin();

		maxval = sci.IDL_MAX(param->vSubArr);

		if (maxval > thresh)
		{
			// SUBTRACT MIN
			for (int k=0;k<xs*ys;++k,++ptr)
			{
				if (*(ptr) < thresh) *(ptr) = 0;
			}

			// COMPUTE CENTROID
			sci.Centroid(param->vSubArr,xs,ys,xc,yc);
		}
		else
		{
			// SET TO REF POSITION
			xc = 0.;
			yc = 0.;
		}

		SH_xtmp[crd] = float(i_xx + xc);
		SH_ytmp[crd] = float(i_yy + yc);
	}

	//QueryPerformanceCounter(&t1);

	//cout << "Execution time: " << double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart) << endl;

//    if ref is True:
//        self.SH_xref = self.SH_xtmp.copy()
//        self.SH_yref = self.SH_ytmp.copy()

//    self.SH_xslp = self.SH_xtmp - self.SH_xref
//    self.SH_yslp = self.SH_ytmp - self.SH_yref

//    self.SH_phot -= self.threshold
//    self.SH_xslp[self.SH_phot < 0] = 0.0
//    self.SH_yslp[self.SH_phot < 0] = 0.0

//    self.SH_comb[0] = self.SH_xslp
//    self.SH_comb[1] = self.SH_yslp
//    self.SH_comb[2] = self.SH_phot / self.SH_phot.max()

//    self.shm_comb.set_data(self.SH_comb)
//    #self.shm_xslp.set_data(self.SH_xslp)
//    #self.shm_yslp.set_data(self.SH_yslp)
//    self.shm_phot_inst.set_data(self.SH_phot)
//    
//    # here is information about the tip-tilt in pixels!
//    # weighted mean version!
//    self.ttx_mean = np.average(self.SH_xslp, weights=self.SH_phot)
//    self.tty_mean = np.average(self.SH_yslp, weights=self.SH_phot)
//    #self.ttx_mean = np.sum(self.SH_xslp * self.SH_phot) / np.sum(self.SH_phot)
//    #self.tty_mean = np.sum(self.SH_yslp * self.SH_phot) / np.sum(self.SH_phot)

//    # original version below
//    #self.ttx_mean = np.median(self.SH_xslp[self.SH_phot > self.threshold])
//    #self.tty_mean = np.median(self.SH_yslp[self.SH_phot > self.threshold])

//    self.update_log()

	return SUCCESS;
}

int _tmain(int argc, _TCHAR* argv[])
{
	MCODE_Export_FITS FITS;

	vector<int> vSize(2);
	vector<int> v_xx;
	vector<int> v_yy;
	vector<int> vImg;
	
	char mychar;

	string data_dir = "C:\\RSI\\MyIDL\\AOC\\simu\\";

	cout << "Starting..." << endl;

	if(FITS.Load(data_dir + "xx.fits",v_xx,vSize) != SUCCESS)
	{
		return SERIOUS;
	}

	if(FITS.Load(data_dir + "yy.fits",v_yy,vSize) != SUCCESS)
	{
		return SERIOUS;
	}

	for (size_t k=0;k<36;++k)
	{
		if(FITS.Load(data_dir + string_format("wfs_cam_%d.fits",k),vImg,vSize) != SUCCESS)
		{
			return SERIOUS;
		}

		LARGE_INTEGER t0;
		LARGE_INTEGER t1;
		LARGE_INTEGER t2;

		//vector<float> SH_xtmp;
		//vector<float> SH_ytmp;
		//calc_SH_data(v_xx,v_yy,vImg,SH_xtmp,SH_ytmp);

		vector<float> SH_xtmp((v_xx.size() - 1)*(v_yy.size() - 1));
		vector<float> SH_ytmp((v_xx.size() - 1)*(v_yy.size() - 1));

		LARGE_INTEGER liFreq;
		QueryPerformanceFrequency(&liFreq);
		QueryPerformanceCounter(&t0);

		calc_SH_data(int(vSize[0]),int(vSize[1]),&vImg[0],int(v_xx.size()),int(v_yy.size()),&v_xx[0],&v_yy[0],&SH_xtmp[0],&SH_ytmp[0]);

		QueryPerformanceCounter(&t1);

		cout << "Execution time: " << double(t1.QuadPart - t0.QuadPart) / double(liFreq.QuadPart) << endl;

		FITS.Save(SH_xtmp,int(SH_xtmp.size()),1,data_dir,string_format("SH_xtmp_%d",k));
		FITS.Save(SH_ytmp,int(SH_ytmp.size()),1,data_dir,string_format("SH_ytmp_%d",k));

		cout << "...Finished." << endl;
	}

	system("pause");

	return 0;
}

