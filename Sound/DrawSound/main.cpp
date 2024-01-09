#include <QApplication>

#include "recordanddrawControl.h"


int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    SoundRecordAndDrawControl control;
    control.show();
    return app.exec();
}
