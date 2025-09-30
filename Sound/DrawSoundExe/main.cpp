//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include <QApplication>

#include "recordanddrawControl.h"
#include "soundmainwindow.h"


int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    SoundMainWindow mw;
    SoundRecordAndDrawControl control(&mw);
    mw.show();
    mw.showMDIChild(&control);
    return app.exec();
}
