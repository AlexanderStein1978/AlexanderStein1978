#include "networkserver.h"
#include "window.h"
#include "logger.h"

#include "QTcpSocket"


NetworkServer::NetworkServer(Window *window, QTcpSocket* socket) : Network(window)
{
    mSocket = socket;
    mTimer.start();
}

void NetworkServer::commandReceived(const Command command)
{
    switch (command)
    {
    case START:
        mWindow->StartCommandReceived();
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
    case RELOAD_POTENTIALS:
        mWindow->emitReloadPotentials();
        break;
    case STOP_CALC:
        mWindow->stopCalc();
        break;
    case ROTATE:
        mWindow->rotate();
        break;
    case DATA_RECEIVED:
        SendData();
        break;
    case GET_SETTINGS_AND_POTENTIALS:
        mWindow->GetSettingsRequestReceived();
        break;
    case ERROR_INCOMPLETE:
    case ERROR_UNKNOWN_COMMAND:
        break;
    default:
        Network::commandReceived(command);
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
    Logger::getLogger().SetNetworkServer(this);
    mTimer.start();
}

void NetworkServer::SendData()
{
    QByteArray data;
    mWindow->copyDataIfNew(data, mCommandMap.key(DATA_FOLLOWING));
    if (data.length() > 28) SendCommand(data);
}

void NetworkServer::SendLogMessage(const QtMsgType type, const QString &time, const QString &function, const QString &file, const QString &message)
{
    static const size_t minSizeLogMessage(44u);
    static const size_t functionOffset(SIZE_OF_COMMAND_STRINGS + 1 + time.size());
    QByteArray string = mCommandMap.key(LOGMESSAGE);
    string.resize(minSizeLogMessage + function.size() + file.size() + message.size());
    string[SIZE_OF_COMMAND_STRINGS] = static_cast<quint8>(type);
    memcpy(string.data() + SIZE_OF_COMMAND_STRINGS + 1, time.data(), time.size());
    size_t offset(functionOffset);
    appendToByteArray(string, offset, function);
    appendToByteArray(string, offset, file);
    appendToByteArray(string, offset, message);
    SendCommand(string);
}

void NetworkServer::appendToByteArray(QByteArray& string, size_t& offset, const QString& stringToAppend)
{
    const quint32 size(stringToAppend.size());
    memcpy(string.data() + offset, &size, sizeof(quint32));
    offset += sizeof(quint32);
    memcpy(string.data() + offset, stringToAppend.data(), size);
    offset += size;
}
