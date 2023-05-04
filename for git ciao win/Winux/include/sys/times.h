#pragma once

//#ifndef _TIMES_H
//
//#include "sys/times.h"
//
//#endif

#ifndef _TIMES_H
#define _TIMES_H

#ifdef COMPILE_MYLIBRARY
  #define MYLIBRARY_EXPORT __declspec(dllexport)
#else
  #define MYLIBRARY_EXPORT __declspec(dllimport)
#endif

#define _WINSOCKAPI_
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock.h>

#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>

MYLIBRARY_EXPORT int gettimeofday(struct timeval* t, void* timezone);

// from linux's sys/times.h

//#include <features.h>

#define __need_clock_t
#include <time.h>

/* Structure describing CPU time used by a process and its children.  */
struct tms
{
    clock_t tms_utime;          /* User CPU time.  */
    clock_t tms_stime;          /* System CPU time.  */

    clock_t tms_cutime;         /* User CPU time of dead children.  */
    clock_t tms_cstime;         /* System CPU time of dead children.  */
};

/* Store the CPU time used by this process and all its
   dead children (and their dead children) in BUFFER.
   Return the elapsed real time, or (clock_t) -1 for errors.
   All times are in CLK_TCKths of a second.  */
MYLIBRARY_EXPORT clock_t times(struct tms* __buffer);

typedef long long suseconds_t;
#endif
#endif