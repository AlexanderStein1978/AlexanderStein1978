#include "logfunction.h"
#include "logger.h"
#include "networkserver.h"

#include <QMessageLogContext>
#include <QDateTime>
#include <QMap>


void LogFunction(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    QString file(QString("%1:%2").arg(context.file ? context.file : "").arg(context.line));
    QString function(context.function ? context.function : "");
    QDateTime dateTime(QDateTime::currentDateTime());
    Logger& logger(Logger::getLogger());
    QString timeString(dateTime.toString("dd.MM.yyyy hh:mm:ss.zzz"));

    NetworkServer* server = logger.GetNetworkServer();
    if (nullptr != server) server->SendLogMessage(type, timeString, function, file, msg);

    logger.LogMessage(QStringList() << timeString << logger.getTypeMap().value(type) << function << file << msg);
}
