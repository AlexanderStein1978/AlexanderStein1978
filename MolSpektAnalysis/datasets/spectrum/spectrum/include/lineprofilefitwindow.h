//
// C++ Interface: LineProfileFitWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINEPROFILEFITWINDOW_H
#define LINEPROFILEFITWINDOW_H


#include "linewindowbase.h"

class QComboBox;
class QLineEdit;
class QPushButton;
class QLabel;

class MainWindow;
class Spektrum;
class Gaussian;


class LineProfileFitWindow : public LineWindowBase
{
    Q_OBJECT
public:
    LineProfileFitWindow(MainWindow *mw, Spektrum* spect = nullptr, Gaussian* line = nullptr);

private slots:
    virtual void LineChanged(const int index) override;
    virtual void SpektrumChanged(const QString &SpectName) override;
    void RunFit();
    void UpdateSigma();

private:
    QLineEdit *MaxIterationEdit, *MinImprovementEdit, *MinFreqEdit, *MaxFreqEdit;
    QPushButton *PerformFitButton;
    QLabel *ResultSigmaLabel;
};

#endif // LINEPROFILEFITWINDOW_H
