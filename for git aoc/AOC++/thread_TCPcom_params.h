#pragma once

#include <cstdlib>
#include <Windows.h>

#include <MCODE.h>
#include <MCODE_Socket.h>

using namespace MultiCODE;

class AOC_Socket : public MCODE_Socket
{
public:
	AOC_Socket();
	virtual ~AOC_Socket();

public:

private:
	virtual void			OnAccept(int nErrorCode);
	virtual void			OnClose(int nErrorCode);
	virtual void			OnConnect(int nErrorCode);
	virtual void			OnSend(int nErrorCode);

};

typedef struct
{
	// THREAD CREATION/TERMINATION EVENTS
	HANDLE			hThreadCreated;
	HANDLE			hTerminateThread;

	// SIGNALING EVENTS
	HANDLE			hDataIncoming;
	HANDLE			hDataProcDone;

	
} ThreadParams_TCPcom;