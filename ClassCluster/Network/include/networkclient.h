#ifndef NETWORKCLIENT_H
#define NETWORKCLIENT_H


#include "network.h"


class NetworkClient : public Network
{
public:
    NetworkClient(Window *window);

    void SendCommand(const Command command, const double parameter);
    void SendCommand(const Command command, const QString parameter);

protected:
    void dataReceived();
};

#endif // NETWORKCLIENT_H
