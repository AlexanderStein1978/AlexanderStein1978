#include <QApplication>

#include "controlwindow.h"
#include "MainWindow.h"


int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    MainWindow MW;
    MW.show();
    MW.showMDIChild(new ControlWindow(&MW));
    app.connect(&app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()));
    return app.exec();
}
