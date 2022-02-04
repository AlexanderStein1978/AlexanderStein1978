#include "include/network.h"


int Network::SIZE_OF_COMMAND_STRINGS = 8;


Network::Network(Window *window) : minimumDataToRead(1), mSocket(nullptr), mWindow(window),
    mCommandMap({"STARTSTA", "STOPSTOP", "RESETRES", "MOVEMOVE", "TRIGSNAP", "WRITSNAP", "RESTSNAP", "SETKINEN", "SETPOTRS", "SETSPEED", "SETSTEPS", "RELOAPOT", "STOPCALC",
                 "ROTATERO", "SETLAYDI", "DATRECED", "DATFOLLO", "ERRORINC", "ERRORUNK"}),
    continueRunning(true)
{
}

Network::~Network()
{
    if (nullptr != mSocket) delete mSocket;
    continueRunning = false;
    wait();
}

void Network::run()
{
    if (nullptr == mSocket) return;
    while (continueRunning)
    {
        if (mSocket->bytesAvailable() >= SIZE_OF_COMMAND_STRINGS) dataReceived();
        msleep(1);
    }
}

void Network::SendCommand(const QByteArray &command)
{
    mSocket->write(command);
    mSocket->flush();
}

void Network::SendCommand(const Command command)
{
    SendCommand(mCommandMap.key(command));
}

