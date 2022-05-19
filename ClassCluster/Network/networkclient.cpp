#include "include/networkclient.h"


#include "window.h"
#include "logger.h"
#include <QTcpSocket>

NetworkClient::NetworkClient(Window *window) : Network(window), mWaitingForData(false), mDataLength(0)
{
    mSocket = new QTcpSocket;
    connect(mSocket, SIGNAL(connected()), this, SIGNAL(IsConnected()));
    connect(mSocket, SIGNAL(disconnected()), this, SIGNAL(ConnectionFailed()));
    connect(mSocket, SIGNAL(error(QAbstractSocket::SocketError)), this, SLOT(ConnectionError()));
}

void NetworkClient::ConnectToServer(const QString IpAddress)
{
    mSocket->connectToHost(IpAddress, 50000);
    mTimer.start();
}

void NetworkClient::ConnectionError()
{
    QAbstractSocket::SocketState state = mSocket->state();
    if (state != QAbstractSocket::ConnectedState && state != QAbstractSocket::HostLookupState && state != QAbstractSocket::ConnectingState) emit ConnectionFailed();
}

void NetworkClient::SendCommand(const Command command, const double parameter)
{
    QByteArray completeCommand;
    const int completeSize = SIZE_OF_COMMAND_STRINGS + sizeof(double);
    completeCommand.reserve(completeSize);
    completeCommand += mCommandMap.key(command);
    qInfo() << "Sending command " << completeCommand << " with value " << parameter;
    completeCommand.resize(completeSize);
    memcpy(completeCommand.data() + SIZE_OF_COMMAND_STRINGS, &parameter, sizeof(double));
    Network::SendCommand(completeCommand);
}

void NetworkClient::SendCommand(const Command command, const QString parameter)
{
    QByteArray data(parameter.toUtf8());
    Network::SendCommand(command, data);
}

void NetworkClient::SendCommand(const Command command, const QByteArray& data)
{
    Network::SendCommand(command, data);
}

void NetworkClient::dataReceived()
{
    if (!mWaitingForData) Network::dataReceived();
    else
    {
        mWaitingForData = false;
        minimumDataToRead = SIZE_OF_COMMAND_STRINGS;
        qint64 dataBytesLength = mDataLength * 24;
        qint64 bytesRead = mSocket->read(mWindow->getDrawingDataToFill(mDataLength), dataBytesLength);
        if (bytesRead < dataBytesLength) Network::SendCommand(ERROR_INCOMPLETE);
        else
        {
            mWindow->update();
            Network::SendCommand(DATA_RECEIVED);
        }
    }
}

void NetworkClient::commandReceived(const Command command)
{
    switch (command)
    {
    case DATA_FOLLOWING:
        {
            bool complete = true;
            double totalEnergy, kineticEnergy(readDouble(complete));
            if (complete) totalEnergy = readDouble(complete);
            if (complete) mDataLength = readUint32(complete);
            if (complete)
            {
                mWindow->updateRemoteEnergies(kineticEnergy, totalEnergy);
                mWaitingForData = true;
                qint64 dataBytesLength = mDataLength * 24;
                if (mSocket->bytesAvailable() >= dataBytesLength) dataReceived();
                else minimumDataToRead = dataBytesLength;
            }
        }
        break;
    case SEND_FLAGS:
        {
            char flags;
            qint64 bytesRead = mSocket->read(&flags, 1);
            if (0u == bytesRead) Network::SendCommand(ERROR_INCOMPLETE);
            else mWindow->flagsReceived(flags);
        }
        break;
    case LOGMESSAGE:
        {
            QtMsgType type;
            if (0u == mSocket->read(reinterpret_cast<char*>(&type), 1))
            {
                Network::SendCommand(ERROR_INCOMPLETE);
                return;
            }
            QString function, file, message;
            QByteArray time;
            if (Network::SIZE_OF_LOG_MESSAGE_TIME > mSocket->read(time.data(), Network::SIZE_OF_LOG_MESSAGE_TIME))
            {
                Network::SendCommand(ERROR_INCOMPLETE);
                return;
            }
            if ((function = readString()).isEmpty()) return;
            if ((file = readString()).isEmpty()) return;
            if ((message = readString()).isEmpty()) return;
            Logger& logger = Logger::getLogger();
            logger.LogRemoteMessage(type, time, function, file, message);
        }
    case ERROR_INCOMPLETE:
    case ERROR_UNKNOWN_COMMAND:
        break;
    default:
        Network::SendCommand(ERROR_UNKNOWN_COMMAND);
        break;
    }
}
