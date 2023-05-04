#pragma once

#include "ImageStruct.h"

#if defined(WINDOWS)
//int clock_gettime(int X, struct timeval *tv);
int myclock_gettime(int, struct mytimespec *spec);
#endif

int ImageCreateSem(IMAGE *image, long NBsem);

int ImageCreate(IMAGE *image, const char *name, long naxis, uint32_t *size, uint8_t atype, int shared, int NBkw);


