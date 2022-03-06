#include "network.h"
#include <QTcpSocket>


const int Network::SIZE_OF_COMMAND_STRINGS = 8;


Network::Network(Window *window) : minimumDataToRead(SIZE_OF_COMMAND_STRINGS), mSocket(nullptr), mWindow(window),
    mCommandMap({{"STARTSTA", START}, {"STOPSTOP", STOP}, {"SEDFLAGS", SEND_FLAGS}, {"RESETRES", RESET}, {"MOVEMOVE", MOVE}, {"TRIGSNAP", TRIGGER_SNAP_SHOT}, {"WRITSNAP", WRITE_SNAP_SHOT},
                 {"RESTSNAP", RESTORE_SNAP_SHOT}, {"SETSETTI", SET_SETTINGS}, {"SETPOTAL", SET_POTENTIAL}, {"SETKINEN", SET_KINETIC_ENERGY}, {"SETPOTRS", SET_POTENTIAL_RANGE_SCALE}, {"SETSPEED", SET_SPEED},
                 {"SETSTEPS", SET_STEP_SIZE}, {"RELOAPOT", RELOAD_POTENTIALS}, {"STOPCALC", STOP_CALC}, {"ROTATERO", ROTATE}, {"SETLAYDI", SET_LAYER_DISTANCE}, {"DATRECED", DATA_RECEIVED},
                 {"DATFOLLO", DATA_FOLLOWING}, {"ERRORINC", ERROR_INCOMPLETE}, {"ERRORUNK", ERROR_UNKNOWN_COMMAND}}),
    continueRunning(true)
{
}

Network::~Network()
{
    continueRunning = false;
    wait();
    if (nullptr != mSocket) delete mSocket;
}

void Network::run()
{
    if (nullptr == mSocket) return;
    while (continueRunning)
    {
        if (mSocket->bytesAvailable() >= minimumDataToRead) dataReceived();
        
        for (auto it = mOldSockets.begin(); it != mOldSockets.end(); ++it)
            if ((*it)->state() == QAbstractSocket::UnconnectedState)
        {
            delete *it;
            *it = nullptr;
            mOldSockets.erase(it);
            break;
        }
        
        msleep(1);
    }
}

void Network::SendCommand(const QByteArray &command)
{
    if (mSocket->state() != QAbstractSocket::ConnectedState) return;
    mSocket->write(command);
    mSocket->flush();
}

void Network::SendCommand(const Command command)
{
    QByteArray string = mCommandMap.key(command);
    qInfo() << "Sending command " << string;
    SendCommand(string);
}

void Network::SendFlags(const char flags)
{
    QByteArray data;
    data.reserve(SIZE_OF_COMMAND_STRINGS + 1);
    data += mCommandMap.key(SEND_FLAGS);
    data[SIZE_OF_COMMAND_STRINGS] = flags;
    qInfo() << "Sending flags: " << static_cast<int>(flags);
    SendCommand(data);
}

void Network::dataReceived()
{
    char buffer[SIZE_OF_COMMAND_STRINGS + 1];
    memset(buffer, 0, SIZE_OF_COMMAND_STRINGS + 1);
    quint64 bytesRead = mSocket->read(buffer, SIZE_OF_COMMAND_STRINGS);
    bool complete = true;
    qInfo() << "Received data: " << buffer;
    if (bytesRead < SIZE_OF_COMMAND_STRINGS)
    {
        SendCommand(ERROR_INCOMPLETE);
        return;
    }
    Command command(mCommandMap.value(buffer, NO_KNOWN_COMMAND));
    if (command == NO_KNOWN_COMMAND)
    {
        SendCommand(ERROR_UNKNOWN_COMMAND);
        RecoverFromError(buffer);
    }
    else commandReceived(command);
}

void Network::RecoverFromError(char *receivedBuffer)
{
    qint64 bytesRead;
    Command command;
    int offset(0);
    char buffer[SIZE_OF_COMMAND_STRINGS + 1];
    memcpy(buffer, receivedBuffer, SIZE_OF_COMMAND_STRINGS + 1);
    do
    {
        for (int i=1; i < SIZE_OF_COMMAND_STRINGS; ++i) buffer[i-1] = buffer[i];
        bytesRead = mSocket->read(buffer + SIZE_OF_COMMAND_STRINGS - 1, 1);
        if (bytesRead == 0) break;
        command = mCommandMap.value(buffer + offset, NO_KNOWN_COMMAND);
    } while (command == NO_KNOWN_COMMAND);
    if (command != NO_KNOWN_COMMAND) commandReceived(command);
}

double Network::readDouble(bool complete)
{
    double value;
    quint64 bytesRead = mSocket->read(reinterpret_cast<char*>(&value), sizeof(double));
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
    quint64 bytesRead = mSocket->read(reinterpret_cast<char*>(&size), sizeof(quint32));
    if (bytesRead < sizeof(quint32))
    {
        SendCommand(ERROR_INCOMPLETE);
        complete = false;
    }
    return size;
}
