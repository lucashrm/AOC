#pragma once

#if defined(WINDOWS)

#ifndef _CRT_NO_TIME_T

	#define _TIMESPEC_DEFINED

//struct timespec
//{
//time_t tv_sec; // Seconds - >= 0
//long tv_nsec; // Nanoseconds - [0, 999999999]
//};

	// WINDOWS IMPLEMENTATION OF clock_gettime
	struct mytimespec
	{
		long tv_sec;
		long tv_nsec;
	};    //header part


#endif

#endif

