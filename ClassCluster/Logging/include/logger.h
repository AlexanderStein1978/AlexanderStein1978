#ifndef LOGGER_H
#define LOGGER_H


#include <QStringList>
#include <QMutex>
#include <QMap>
#include <QObject>
#include "loglist.h"


class LogWindow;


class Logger : public QObject
{
    Q_OBJECT
public:
    static Logger& getLogger();

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void SetLogWindow(LogWindow* window);
    void LogMessage(QStringList& message);
    void LogRemoteMessage(const QtMsgType type, const QString& time, const QString& function, const QString& file, const QString& message);

    inline void SendLogMessageToClient(const QtMsgType type, const QString& time, const QString& function, const QString& file, const QString& message)
    {
        emit SendLogMessage(type, time, function, file, message);
    }

    inline const QMap<QtMsgType, QString>& getTypeMap()
    {
        return mTypeMap;
    }

signals:
    void SendLogMessage(const QtMsgType type, const QString& time, const QString& function, const QString& file, const QString& message);

private:
    Logger();
    virtual ~Logger();

    const QMap<QtMsgType, QString> mTypeMap;
    LogList mMessageBuffer;
    QMutex mMutex;
    LogWindow* mLogWindow;
};

#endif // LOGGER_H
