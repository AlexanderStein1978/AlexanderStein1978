#ifndef NETWORKCLIENT_H
#define NETWORKCLIENT_H


#include "network.h"


class NetworkClient : public Network
{
    Q_OBJECT

public:
    NetworkClient(Window *window);

    void ConnectToServer(const QString IpAddress);
    void SendCommand(const Command command, const double parameter);
    void SendCommand(const Command command, const QString parameter);
    void SendCommand(const Command command, const QByteArray& data);

    inline void SendCommand(const Command command)
    {
        Network::SendCommand(command);
    }

signals:
    void IsConnected();
    void ConnectionFailed();

protected:
    void dataReceived() override;
    void commandReceived(const Command command) override;

private slots:
    void ConnectionError();

private:
    void readData();

    bool mWaitingForData;
    qint32 mDataLength;
};

#endif // NETWORKCLIENT_H
