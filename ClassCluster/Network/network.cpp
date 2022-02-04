#include "include/network.h"
#include <QTcpSocket>


int Network::SIZE_OF_COMMAND_STRINGS = 8;


Network::Network(Window *window) : minimumDataToRead(SIZE_OF_COMMAND_STRINGS), mSocket(nullptr), mWindow(window),
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
        if (mSocket->bytesAvailable() >= minimumDataToRead) dataReceived();
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

void Network::dataReceived()
{
    char buffer[SIZE_OF_COMMAND_STRINGS + 1];
    memset(buffer, 0, SIZE_OF_COMMAND_STRINGS + 1);
    quint64 bytesRead = mSocket->readData(buffer, SIZE_OF_COMMAND_STRINGS);
    bool complete = true;
    if (bytesRead < SIZE_OF_COMMAND_STRINGS)
    {
        SendCommand(ERROR_INCOMPLETE);
        return;
    }
    Command command(mCommandMap.value(buffer, ERROR_UNKNOWN_COMMAND));
    commandReceived(command);
}

double Network::readDouble(bool complete)
{
    double value;
    bytesRead = mSocket->readData(reinterpret_cast<char*>(&value), sizeof(double));
    if (bytesRead < sizeof(double))
    {
        SendCommand(ERROR_INCOMPLETE);
        complete = false;
    }
    return value;
}

quint32 Network::readUint32(bool complete)
{
    quint32 size;
    bytesRead = mSocket->readData(reinterpret_cast<char*>(&size), sizeof(quint32));
    if (bytesRead < sizeof(quint32))
    {
        SendCommand(ERROR_INCOMPLETE);
        complete = false;
    }
    return size;
}
