//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include <QApplication>

#include "controlwindow.h"
#include "MainWindow.h"
#include "logfunction.h"
#include "logwindow.h"
#include "logger.h"


int main(int argc, char **argv)
{
    qInstallMessageHandler(LogFunction);
    QApplication app(argc, argv);
    MainWindow MW;
    MW.show();
    MW.showMDIChild(new ControlWindow(&MW));
    LogWindow* logwindow = new LogWindow(&MW);
    MW.showMDIChild(logwindow);
    Logger::getLogger().SetLogWindow(logwindow);
    app.connect(&app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()));
    return app.exec();
}
