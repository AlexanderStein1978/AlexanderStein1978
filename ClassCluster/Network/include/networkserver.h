#ifndef NETWORKSERVER_H
#define NETWORKSERVER_H


#include "network.h"


class NetworkServer : public Network
{
public:
    NetworkServer(Window *window, QTcpSocket *socket);

    bool IsConnected();
    void NewConnection(QTcpSocket *socket);
    void SendData();

protected:
    void commandReceived(const Command command) override;

private:
    QString readString();
};

#endif // NETWORKSERVER_H
