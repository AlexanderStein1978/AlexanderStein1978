#include "logger.h"
#include "logwindow.h"
#include "networkserver.h"
#include <QMutexLocker>


Logger::Logger() 
: mTypeMap({{QtDebugMsg, "Debug"}, {QtInfoMsg, "Info"}, {QtWarningMsg, "Warning"}, {QtCriticalMsg, "Critical"}, {QtFatalMsg, "Fatal"}})
, mLogWindow(nullptr), mNetworkServer(nullptr)
{
}

Logger& Logger::getLogger()
{
    static Logger instance;
    return instance;
}

void Logger::LogMessage(QStringList &message)
{
    QMutexLocker lock(&mMutex);
    if (nullptr != mLogWindow) mLogWindow->LogMessage(message);
    else mMessageBuffer << message;
}

void Logger::LogRemoteMessage(const QtMsgType type, const QString &time, const QString &function, const QString &file, const QString &message)
{
    QStringList finalMessage;
    finalMessage << time << ("REMOTE " + mTypeMap.value(type)) << function << file << message;
    LogMessage(finalMessage);
}

void Logger::SetLogWindow(LogWindow *window)
{
    QMutexLocker lock(&mMutex);
    mLogWindow = window;
    window->SetMessageBuffer(mMutex, mMessageBuffer);
}

NetworkServer* Logger::GetNetworkServer() const
{
    if (nullptr != mNetworkServer && mNetworkServer->IsConnected()) return mNetworkServer;
    else return nullptr;
}
