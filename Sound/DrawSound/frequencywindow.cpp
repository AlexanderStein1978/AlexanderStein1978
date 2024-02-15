#include "frequencywindow.h"
#include "fit.h"
#include "utils.h"
#include "datensatz.h"
#include "soundwindow.h"
#include "windowselectdialog.h"
#include "recordanddrawControl.h"

#include <QAction>
#include <QMenu>


FrequencyWindow::FrequencyWindow(SoundRecordAndDrawControl *const control, const int sampleRate) : SoundDrawWindow(control, sampleRate, 0)
{
    mRescaleOnSetAndAdd = false;
    QAction *backTransformAct = new QAction("Back transform", this), *copyToWindowAct = new QAction("Copy to window...", this);
    mPopupMenu->addAction(backTransformAct);
    mPopupMenu->addAction(copyToWindowAct);
    connect(backTransformAct, SIGNAL(triggered()), this, SLOT(BackTransform()));
    connect(copyToWindowAct, SIGNAL(toggled(bool)), this, SLOT(CopyAllDataToNewWindow(bool)));
}

void FrequencyWindow::BackTransform()
{
    int xStart, xStop, n, i, inputLength = getSoundDataRange(xStart, xStop), transLength = inputLength + 1;
    double *realInputdata = new double[inputLength], *imaginaryInputdata = new double[inputLength], **realTransData = Create(transLength, 2), **imaginaryTransData = Create(transLength, 2), step(0.0);
    QString FName("Backtransformed.dat");
    for (n = xStart, i=0; n <= xStop; ++n, ++i)
    {
        realInputdata[i] = Daten[0].GetValue(n, 1);
        imaginaryInputdata[i] = Daten[1].GetValue(n, 1);
        if (step == 0.0 && n>0) step = Daten[0].GetValue(n, 0) - Daten[0].GetValue(n-1, 0);
    }
    backtransformFFT(realInputdata, imaginaryInputdata, inputLength, step, realTransData, imaginaryTransData);
    SoundWindow* window = new SoundWindow(mControl, FName, mSampleRate);
    window->setData(realTransData, transLength);
    window->addData(imaginaryTransData, transLength);
    window->show();
    delete[] realInputdata;
    delete[] imaginaryInputdata;
    Destroy(realTransData, transLength);
    Destroy(imaginaryTransData, transLength);
}

void FrequencyWindow::CopyAllDataToNewWindow()
{
    WindowSelectDialog windowDialog(mControl->GetTheFrequencyWindows(), this);
    if (windowDialog.exec() == QDialog::Rejected) return;
    int length = Daten->GetDSL(), index = windowDialog.GetSelection();
    DiagWindow* window;
    if (0 > index && index >= mControl->GetNumberFrequencyWindows())
    {
        NameSelectionDialog nameDialog(this);
        if (nameDialog.exec() == QDialog::Rejected) return;
        window = new DiagWindow;
        window->setUnits("Frequency [Hz]", "Intensity");
        window->setName(nameDialog.GetName());
        mControl->AddFrequencyWindow(window);
    }
    else window = mControl->GetWindow(index);
    double **realData = Create(length, 2), **imaginaryData = Create(length, 2);
    for (int n=0; n < length; ++n)
    {
        realData[n][0] = Daten[0].GetValue(n, 0);
        realData[n][1] = Daten[0].GetValue(n, 1);
        imaginaryData[n][0] = Daten[1].GetValue(n, 0);
        imaginaryData[n][1] = Daten[1].GetValue(n, 1);
    }
    window->addData(realData, length);
    window->addData(imaginaryData, length);
    if (!window->isVisible()) window->show();
    Destroy(realData, length);
    Destroy(imaginaryData, length);
}
