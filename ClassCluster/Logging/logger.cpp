#include "logger.h"
#include "logwindow.h"
#include <QMutexLocker>


Logger::Logger() : mLogWindow(nullptr)
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

void Logger::SetLogWindow(LogWindow *window)
{
    QMutexLocker lock(&mMutex);
    mLogWindow = window;
    window->SetMessageBuffer(mMutex, mMessageBuffer);
}
