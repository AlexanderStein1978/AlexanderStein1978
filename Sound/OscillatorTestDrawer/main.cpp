//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include <QApplication>

#include "oscillator.h"
#include "DiagWindow.h"
#include "utils.h"


int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    Oscillator oscillator;
    oscillator.initialize(20, 2, 1e-2);
    double t, **data = Create(3000, 2);
    int n;
    oscillator.newValue(-10.0);
    for (n=0, t=0.0; n < 3000; ++n, t += 1e-2)
    {
        data[n][0] = t;
        oscillator.newValue(0.0);
        data[n][1] = oscillator.getLastAmp(); //      (n < 100 ? oscillator.getAmplitudeFor(t) : oscillator.getAmplitudeFor(t - 1.0));
    }
    DiagWindow window;
    window.setData(data, 3000);
    window.show();
    return app.exec();
}
