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
    void dataReceived() override;
    void commandReceived(const Command command) override;

private:
    void readData();

    bool mWaitingForData;
    qint32 mDataLength;
};

#endif // NETWORKCLIENT_H
