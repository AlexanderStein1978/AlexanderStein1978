//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

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
    QStringList message;
    message << timeString << logger.getTypeMap().value(type) << function << file << msg;
    logger.LogMessage(message);
    logger.SendLogMessageToClient(type, timeString, function, file, msg);
}
