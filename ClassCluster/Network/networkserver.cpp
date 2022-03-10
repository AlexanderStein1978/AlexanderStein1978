#include "networkserver.h"
#include "window.h"
#include "QTcpSocket"


NetworkServer::NetworkServer(Window *window, QTcpSocket* socket) : Network(window)
{
    mSocket = socket;
    mTimer.start();
}

void NetworkServer::commandReceived(const Command command)
{
    bool complete(true);
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

bool NetworkServer::IsConnected()
{
    return (mSocket->state() == QAbstractSocket::ConnectedState);
}

void NetworkServer::NewConnection(QTcpSocket *socket)
{
    mTimer.stop();
    mSocket->disconnectFromHost();
    mOldSockets.push_back(mSocket);
    mSocket = socket;
    mTimer.start();
}

QString NetworkServer::readString()
{
    bool complete = true;
    quint32 size(readUint32(complete));
    if (!complete) return "";
    QByteArray buffer;
    buffer.resize(size);
    qint64 bytesRead = mSocket->read(buffer.data(), size);
    if (bytesRead < size)
    {
        SendCommand(ERROR_INCOMPLETE);
        return "";
    }
    return buffer;
}

void NetworkServer::SendData()
{
    QByteArray data;
    mWindow->copyDataIfNew(data, mCommandMap.key(DATA_FOLLOWING));
    if (data.length() > 28) SendCommand(data);
}
