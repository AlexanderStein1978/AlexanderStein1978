#ifndef LOGGER_H
#define LOGGER_H


#include <QStringList>
#include <QMutex>
#include <QMap>


class LogWindow;
class NetworkServer;


class Logger
{
public:
    static Logger& getLogger();

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void SetLogWindow(LogWindow* window);
    void LogMessage(QStringList& message);
    void LogRemoteMessage(const QtMsgType type, const QString& time, const QString& function, const QString& file, const QString& message);
    NetworkServer* GetNetworkServer() const;

    inline void SetNetworkServer(NetworkServer* server)
    {
        mNetworkServer = server;
    }

    inline const QMap<QtMsgType, QString>& getTypeMap()
    {
        return mTypeMap;
    }

private:
    Logger();

    const QMap<QtMsgType, QString> mTypeMap;
    QList<QStringList> mMessageBuffer;
    QMutex mMutex;
    LogWindow* mLogWindow;
    NetworkServer* mNetworkServer;
};

#endif // LOGGER_H
