#include "include/networkclient.h"


#include "window.h"
#include <QTcpSocket>

NetworkClient::NetworkClient(Window *window) : Network(window), mWaitingForData(false), mDataLength(0)
{
}

void NetworkClient::SendCommand(const Command command, const double parameter)
{
    QByteArray completeCommand;
    const int completeSize = SIZE_OF_COMMAND_STRINGS + sizeof(double);
    completeCommand.reserve(completeSize);
    completeCommand += mCommandStrings[command];
    completeCommand.resize(completeSize);
    memcpy(completeCommand.data() + command.size(), &parameter, sizeof(double));
    SendCommand(completeCommand);
}

void NetworkClient::SendCommand(const Command command, const QString parameter)
{
    QByteArray completeCommand;
    quint32 size = parameter.size();
    int fullSize = SIZE_OF_COMMAND_STRINGS + sizeof(quint32) + size;
    completeCommand.reserve(fullSize);
    completeCommand += mCommandStrings[command];
    completeCommand.resize(command.size() + sizeof(quint32));
    memcpy(completeCommand.data() + command.size(), &size, sizeof(quint32));
    completeCommand += parameter;
    SendCommand(completeCommand);
}

void NetworkClient::dataReceived()
{
    if (!mWaitingForData) Network::dataReceived();
    else
    {
        mWaitingForData = false;
        minimumDataToRead = SIZE_OF_COMMAND_STRINGS;
        qint64 dataBytesLength = mDataLength * 24;
        qint64 bytesRead = mSocket->readData(mWindow->getDrawingDataToFill(mDataLength), dataBytesLength);
        if (bytesRead < dataBytesLength) SendCommand(ERROR_INCOMPLETE);
        else
        {
            mWindow->update();
            SendCommand(DATA_RECEIVED);
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
                if (mSocket->bytesAvailable() >= dataLength) dataReceived();
                else minimumDataToRead = dataBytesLength;
            }
        }
        break;
    case ERROR_INCOMPLETE:
    case ERROR_UNKNOWN_COMMAND:
        break;
    default:
        SendCommand(ERROR_UNKNOWN_COMMAND);
        break;
    }
}
