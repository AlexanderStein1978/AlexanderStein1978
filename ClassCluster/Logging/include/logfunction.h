//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LOGFUNCTION_H
#define LOGFUNCTION_H


#include <QMessageLogContext>


void LogFunction(QtMsgType type, const QMessageLogContext &context, const QString &msg);


#endif // LOGFUNCTION_H
