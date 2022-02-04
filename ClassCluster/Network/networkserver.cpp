#include "networkserver.h"
#include "window.h"


NetworkServer::NetworkServer(Window *window) : Network(window)
{
}

void NetworkServer::dataReceived()
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
    switch (command)
    {
    case START:
        mWindow->start();
        SendData();
        break;
    case STOP:
        mWindow->stop();
        break;
    case RESET:
        mWindow->reset();
        break;
    case MOVE:
        mWindow->move();
        break;
    case TRIGGER_SNAP_SHOT:
        mWindow->triggerSnapShot();
        break;
    case WRITE_SNAP_SHOT:
        {
            QString filename = readString();
            if (!filename.isEmpty()) mWindow->writeSnapShot(filename);
        }
        break;
    case RESTORE_SNAP_SHOT:
        {
            QString filename = readString();
            bool isMoving;
            if (!filename.isEmpty()) mWindow->restoreSnapShot(isMoving, filename);
        }
        break;
    case SET_KINETIC_ENERGY:
        {
            double E = readDouble(complete);
            if (complete) mWindow->setKineticEnergy(E);
        }
        break;
    case SET_POTENTIAL_RANGE_SCALE:
        {
            double Scale = readDouble(complete);
            if (complete) mWindow->setPotentialRangeScale(Scale);
        }
        break;
    case SET_SPEED:
        {
            double Speed = readDouble(complete);
            if (complete) mWindow->setSpeed(Speed);
        }
        break;
    case SET_STEP_SIZE:
        {
            double Size = readDouble(complete);
            if (complete) mWindow->setStepSize(Size);
        }
        break;
    case RELOAD_POTENTIALS:
        mWindow->emitReloadPotentials();
        break;
    case STOP_CALC:
        mWindow->stopCalc();
        break;
    case ROTATE:
        mWindow->rotate();
        break;
    case SET_LAYER_DISTANCE:
        {
            double distance = readDouble(complete);
            if (complete) mWindow->setStepSize(distance);
        }
        break;
    case DATA_RECEIVED:
        SendData();
        break;
    case ERROR_INCOMPLETE:
    case ERROR_UNKNOWN_COMMAND:
        break;
    default:
        SendCommand(ERROR_UNKNOWN_COMMAND);
        break;
    }
}

QString NetworkServer::readString()
{
    quint32 size;
    bytesRead = mSocket->readData(reinterpret_cast<char*>(&size), sizeof(quint32));
    if (bytesRead < quint32)
    {
        SendCommand(ERROR_INCOMPLETE);
        return "";
    }
    QString filename;
    filename.resize(size);
    bytesRead = mSocket->readData(filename.data(), size);
    if (bytesRead < size)
    {
        SendCommand(ERROR_INCOMPLETE);
        return "";
    }
    return filename;
}

double NetworkServer::readDouble(bool complete)
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

void NetworkServer::SendData()
{
    QByteArray data;
    mWindow->copyDataIfNew(data, mCommandMap.key(DATA_FOLLOWING));
    if (data.length() > 24) SendCommand(data);
}
