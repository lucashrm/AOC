#include <stdio.h>
#include <getpagesize.h>
#include <sys/times.h>

// Code for determining page size swiped from Python's mmapmodule.c
#if defined(HAVE_SYSCONF) && defined(_SC_PAGESIZE)
static int
my_getpagesize(void)
{
	return sysconf(_SC_PAGESIZE);
}
#else
#include <unistd.h>
#define my_getpagesize getpagesize
#endif

static int getpagesize(void) {
    static int cachedPageSize = 0;
    if (cachedPageSize == 0) {
        SYSTEM_INFO si;
        GetSystemInfo(&si);
        cachedPageSize = si.dwPageSize;
    }
    return cachedPageSize;
}

int main(void) { 
    printf("%d\n", getpagesize());

    //struct tms tmy;

    
    return 0; 
}
