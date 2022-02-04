#include "include/networkclient.h"

#include <QTcpSocket>

NetworkClient::NetworkClient(Window *window) : Network(window)
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
