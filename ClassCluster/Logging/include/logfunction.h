#ifndef LOGFUNCTION_H
#define LOGFUNCTION_H


#include <QMessageLogContext>


void LogFunction(QtMsgType type, const QMessageLogContext &context, const QString &msg);


#endif // LOGFUNCTION_H
