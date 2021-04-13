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
class QCheckBox;

class MainWindow;
class Spektrum;
class Gaussian;
class LineDialog;


class LineProfileFitWindow : public LineWindowBase
{
    Q_OBJECT
public:
    LineProfileFitWindow(MainWindow *mw, Spektrum* spect = nullptr, LineProfile* line = nullptr);

    inline void setLineDialog(LineDialog* const Dialog)
    {
        mLineDialog = Dialog;
    }

private slots:
    void RunFit();
    void UpdateSigma();
    void RangeEdited();
    void RangeChanged(const double newMinE, const double newMaxE);

private:
    void disconnectSpectrum() override;
    void connectSpectrum() override;
    void lineChanged() override;
    
    QLineEdit *MaxIterationEdit, *MinImprovementEdit, *MinFreqEdit, *MaxFreqEdit;
    QComboBox *LineTypeBox;
    QPushButton *PerformFitButton;
    QLabel *ResultSigmaLabel;
    LineDialog* mLineDialog;
    QCheckBox* mWithSaturationBox;
};

#endif // LINEPROFILEFITWINDOW_H
