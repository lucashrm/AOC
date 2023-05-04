/*
    Copyright (c) 2011, Dongsheng Song <songdongsheng@live.cn>
	
    Licensed to the Apache Software Foundation (ASF) under one or more
    contributor license agreements.  See the NOTICE file distributed with
    this work for additional information regarding copyright ownership.
    The ASF licenses this file to You under the Apache License, Version 2.0
    (the "License"); you may not use this file except in compliance with
    the License.  You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/

/*
	Simple Windows replacement for POSIX semaphores 
	Modified by Daniel Tillett from libpthread <http://github.com/songdongsheng/libpthread>
	Copyright (c) 2015, Daniel Tillett <daniel.tillett @ gmail.com>
*/
#ifdef COMPILE_MYLIBRARY
#define MYLIBRARY_EXPORT __declspec(dllexport)
#else
#define MYLIBRARY_EXPORT __declspec(dllimport)
#endif


#ifndef _SEMAPHORE_H_
#define _SEMAPHORE_H_   1

/**
    @file semaphore.h
    @brief POSIX Semaphore Definitions and Routines
*/

/**
    @defgroup sem POSIX Semaphore Definitions and Routines
    @{
*/

#include <errno.h> /* Adding definition of EINVAL, ETIMEDOUT, ..., etc. */
#include <fcntl.h> /* Adding O_CREAT definition. */
#include <stdio.h>
#define WIN32_MEAN_AND_LEAN
#include <Windows.h>
#include <winsock.h>

#ifndef PTHREAD_PROCESS_SHARED
#define PTHREAD_PROCESS_PRIVATE	0
#define PTHREAD_PROCESS_SHARED	1
#endif

/* Support POSIX.1b semaphores.  */
#ifndef _POSIX_SEMAPHORES
#define _POSIX_SEMAPHORES       200809L
#endif

#ifndef SEM_VALUE_MAX
#define SEM_VALUE_MAX           INT_MAX
#endif

#ifndef MY_SEM_FAILED
#define MY_SEM_FAILED              NULL
#endif

#define UNUSED(x)				(void)(x)

#ifndef ETIMEDOUT
#define ETIMEDOUT				138 /* This is the value in VC 2010. */
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef void  *my_sem_t;

#ifndef _MODE_T_
typedef unsigned short _mode_t;
#define _MODE_T_  1

#ifndef NO_OLDNAMES
typedef _mode_t mode_t;
#endif
#endif  /* _MODE_T_ */

typedef struct {
	HANDLE handle;
	} arch_sem_t;

//#if !defined(_TIMESPEC_DEFINED)
//struct t
// imespec {
//	time_t  tv_sec;       /* Seconds */
//	long    tv_nsec;      /* Nanoseconds */
//	};
//
//struct itimerspec {
//	struct timespec  it_interval; /* Timer period */
//	struct timespec  it_value;    /* Timer expiration */
//	};
//#define _TIMESPEC_DEFINED       1
//#else
//#include <time.h>
//#endif  /* _TIMESPEC_DEFINED */

MYLIBRARY_EXPORT int my_sem_init(my_sem_t *sem, int pshared, unsigned int value);
MYLIBRARY_EXPORT int my_sem_wait(my_sem_t*sem);
MYLIBRARY_EXPORT int my_sem_trywait(my_sem_t*sem);
MYLIBRARY_EXPORT int my_sem_timedwait(my_sem_t*sem, const struct timespec *abs_timeout);
MYLIBRARY_EXPORT int my_sem_post(my_sem_t*sem);
MYLIBRARY_EXPORT int my_sem_getvalue(my_sem_t*sem, int *value);
MYLIBRARY_EXPORT int my_sem_destroy(my_sem_t*sem);
MYLIBRARY_EXPORT my_sem_t*my_sem_open(const char *name, int oflag, mode_t mode, unsigned int value);
MYLIBRARY_EXPORT int my_sem_close(my_sem_t*sem);
MYLIBRARY_EXPORT int my_sem_unlink(const char *name);

#ifdef __cplusplus
	}
#endif

/** @} */

#endif /* _SEMAPHORE_H_ */
