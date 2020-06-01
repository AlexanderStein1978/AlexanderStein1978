#include <QApplication>

#include "controlwindow.h"

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    ControlWindow w;
    w.show();
	app.connect(&app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()));
    return app.exec();
}
