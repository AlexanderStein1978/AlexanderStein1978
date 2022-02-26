#ifndef LOGGER_H
#define LOGGER_H


#include <QStringList>
#include <QMutex>


class LogWindow;


class Logger
{
public:
    static Logger& getLogger();

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void SetLogWindow(LogWindow* window);
    void LogMessage(QStringList& message);

private:
    Logger();

    QList<QStringList> mMessageBuffer;
    QMutex mMutex;
    LogWindow* mLogWindow;
};

#endif // LOGGER_H
