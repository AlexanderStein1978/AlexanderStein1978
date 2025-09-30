//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef POTENTIALDEFINER_H
#define POTENTIALDEFINER_H


#include "DiagWindow.h"
#include "potentialdefinerinputdata.h"
#include "vector.h"


class Window;


class PotentialDefiner : public DiagWindow
{
    Q_OBJECT
public:
    PotentialDefiner(Window* window, MainWindow *MW);

protected:
    void PSpektrum(QPainter &P, const QRect & A, bool PrintFN) override;

private slots:
    void Calculate();

private:
    PotentialDefinerInputData mData;
    Window* mWindow;
    Vector mStart, mEnd;
    QPushButton *mCalcB;
    QLineEdit *mParticleIndexInput, *mDirectionXInput, *mDirectionYInput, *mDirectionZInput;
};

#endif // POTENTIALDEFINER_H
