#include "logfunction.h"
#include "logger.h"

#include <QMessageLogContext>
#include <QDateTime>
#include <QMap>


void LogFunction(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    static QMap<QtMsgType, QString> typeMap({{QtDebugMsg, "Debug"}, {QtInfoMsg, "Info"}, {QtWarningMsg, "Warning"}, {QtCriticalMsg, "Critical"}, {QtFatalMsg, "Fatal"}});
    QString file(QString("%s:%d").arg(context.file ? context.file : "").arg(context.line));
    QString function(context.function ? context.function : "");
    QDateTime dateTime(QDateTime::currentDateTime());
    Logger& logger(Logger::getLogger());

    logger.LogMessage(QStringList() << dateTime.toString("dd.MM.yyyy hh:mm:ss.zzz") << typeMap.value(type) << function << file << msg);
}
