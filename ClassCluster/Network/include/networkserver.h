#ifndef NETWORKSERVER_H
#define NETWORKSERVER_H


#include "network.h"


class NetworkServer : public Network
{
public:
    NetworkServer(Window *window);

    void SendData();

protected:
    void dataReceived();

private:
    QString readString();
    double readDouble(bool complete);
};

#endif // NETWORKSERVER_H
