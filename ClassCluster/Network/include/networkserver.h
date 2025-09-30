//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef NETWORKSERVER_H
#define NETWORKSERVER_H


#include "network.h"


class NetworkServer : public Network
{
    Q_OBJECT
public:
    NetworkServer(Window *window, QTcpSocket *socket);

    bool IsConnected();
    void NewConnection(QTcpSocket *socket);
    void SendData();

private slots:
    void SendLogMessage(const QtMsgType type, const QString& time, const QString& function, const QString& file, const QString& message);

protected:
    void commandReceived(const Command command) override;
    void appendToByteArray(QByteArray& string, size_t& offset, const QString& stringToAppend);
};

#endif // NETWORKSERVER_H
