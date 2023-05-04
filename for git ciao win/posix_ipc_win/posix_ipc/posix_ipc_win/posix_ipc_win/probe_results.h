/*
This header file was generated when you ran setup. Once created, the setup
process won't overwrite it, so you can adjust the values by hand and
recompile if you need to.

On your platform, this file may contain only this comment -- that's OK!

To enable lots of debug output, add this line and re-run setup.py:
#define POSIX_IPC_DEBUG

To recreate this file, just delete it and re-run setup.py.
*/

#define POSIX_IPC_VERSION		"1.0.4"
#define SEM_GETVALUE_EXISTS		
#define SEM_TIMEDWAIT_EXISTS		
#define MESSAGE_QUEUE_SUPPORT_EXISTS		
#define QUEUE_MESSAGES_MAX_DEFAULT		1024
#define QUEUE_MESSAGE_SIZE_MAX_DEFAULT		8192
#define QUEUE_PRIORITY_MAX		31U
#ifndef PAGE_SIZE
#define PAGE_SIZE		4096
#endif
