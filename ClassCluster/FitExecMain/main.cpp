#include <QApplication>

#include "fitexeccontrol.h"


int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    FitExecControl control;
    return app.exec();
}
