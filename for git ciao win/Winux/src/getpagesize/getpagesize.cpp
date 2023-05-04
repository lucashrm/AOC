#include <Windows.h>
#include <getpagesize.h>

int getpagesize(void) {
  static int cachedPageSize = 0;
  if (cachedPageSize == 0) {
      SYSTEM_INFO si;
      GetSystemInfo(&si);
      cachedPageSize =  si.dwPageSize;
  }
  return cachedPageSize;
}