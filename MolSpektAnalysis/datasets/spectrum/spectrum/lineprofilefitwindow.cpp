//
// C++ Implementation: LineProfileFitWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "lineprofilefitwindow.h"
#include "MainWindow.h"
#include "Spektrum.h"
#include "gaussian.h"
#include "linedialog.h"

#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QGridLayout>
#include <QIntValidator>
#include <QDoubleValidator>
#include <QCheckBox>


LineProfileFitWindow::LineProfileFitWindow(MainWindow* mw, Spektrum *spect, Gaussian *line) : LineWindowBase(mw, spect, line),
    MaxIterationEdit(new QLineEdit("100", this)), MinImprovementEdit(new QLineEdit("0.01", this)), MinFreqEdit(new QLineEdit(this)), MaxFreqEdit(new QLineEdit(this)),
    PerformFitButton(new QPushButton("Run fit", this)), ResultSigmaLabel(new QLabel(this)), mLineDialog(nullptr), mWithSaturationBox(new QCheckBox("Consider saturation", this))
{
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Spectrum:", this), 0, 0);
    L->addWidget(SpektrumBox, 0, 1);
    L->addWidget(new QLabel("Line:", this), 0, 2);
    L->addWidget(LineBox, 0, 3);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Max iterations:", this), 2, 0);
    L->addWidget(MaxIterationEdit, 2, 1);
    L->addWidget(new QLabel("Min improvement:", this), 2, 2);
    L->addWidget(MinImprovementEdit, 2, 3);
    L->addWidget(new QLabel("Min frequency:", this), 3, 0);
    L->addWidget(MinFreqEdit, 3, 1);
    L->addWidget(new QLabel("Max frequency:", this), 3, 2);
    L->addWidget(MaxFreqEdit, 3, 3);
    L->setRowMinimumHeight(4, 20);
    L->addWidget(mWithSaturationBox, 5, 0, 1, 2);
    L->addWidget(PerformFitButton, 5, 2, 1, 2);
    L->addWidget(ResultSigmaLabel, 6, 0, 1, 4);
    MaxIterationEdit->setValidator(new QIntValidator(1, 1000000, MaxIterationEdit));
    MinImprovementEdit->setValidator(new QDoubleValidator(1e-10, 1e10, 0, MinImprovementEdit));
    MinFreqEdit->setValidator(new QDoubleValidator(0.001, 1e10, 0, MinFreqEdit));
    MaxFreqEdit->setValidator(new QDoubleValidator(0.001, 1e10, 0, MaxFreqEdit));
    connect(SpektrumBox, SIGNAL(currentTextChanged(QString)), this, SLOT(SpektrumChanged(QString)));
    connect(LineBox, SIGNAL(currentIndexChanged(int)), this, SLOT(LineChanged(int)));
    connect(MinFreqEdit, SIGNAL(editingFinished()), this, SLOT(RangeEdited()));
    connect(MaxFreqEdit, SIGNAL(editingFinished()), this, SLOT(RangeEdited()));
    connect(PerformFitButton, SIGNAL(clicked()), this, SLOT(RunFit()));
    focusInEvent(nullptr);
}

void LineProfileFitWindow::lineChanged()
{
    UpdateSigma();
}

void LineProfileFitWindow::RunFit()
{
    if (nullptr != mSpektrum)
    {
        int lineIndex = LineBox->currentIndex();
        ResultSigmaLabel->setText("Sigma: " + QString::number(
            mSpektrum->FitGaussianLineProfile(lineIndex, mWithSaturationBox->isChecked(), MaxIterationEdit->text().toDouble(), MinImprovementEdit->text().toDouble(),
                                              MinFreqEdit->text().toDouble(), MaxFreqEdit->text().toDouble())));
        if (nullptr == mLineDialog)
        {
            mLineDialog = new LineDialog(MW, mSpektrum, mSpektrum->GetFittedLine(lineIndex));
            MW->showMDIChild(mLineDialog);
        }
        else mLineDialog->selectLine(lineIndex);
    }
}

void LineProfileFitWindow::connectSpectrum()
{
    if (nullptr != mSpektrum)
    {
        LineBox->addItem("new");
        connect(mSpektrum, SIGNAL(NumberOfFittedLinesChanged()), this, SLOT(NumberOfLinesChanged()));
        connect(mSpektrum, SIGNAL(propertiesChanged()), this, SLOT(UpdateSigma()));
        connect(mSpektrum, SIGNAL(SelectedRangeChanged(double,double)), this, SLOT(RangeChanged(double,double)));
    }
}

void LineProfileFitWindow::disconnectSpectrum()
{
    if (nullptr != mSpektrum)
    {
        disconnect(mSpektrum, SIGNAL(NumberOfFittedLinesChanged()), this, SLOT(NumberOfLinesChanged()));
        disconnect(mSpektrum, SIGNAL(propertiesChanged()), this, SLOT(UpdateSigma()));
        disconnect(mSpektrum, SIGNAL(SelectedRangeChanged(double,double)), this, SLOT(RangeChanged(double,double)));
    }
}

void LineProfileFitWindow::UpdateSigma()
{
    if (nullptr != mLine)
    {
        double EStart, EEnd;
        mLine->GetDataRange(EStart, EEnd);
        MinFreqEdit->setText(QString::number(EStart, 'f', 6));
        MaxFreqEdit->setText(QString::number(EEnd, 'f', 6));
        ResultSigmaLabel->setText("Sigma: " + QString::number(mLine->GetSigma()));
        mWithSaturationBox->setChecked(mLine->isWithSaturation());
        mWithSaturationBox->setEnabled(false);
    }
    else mWithSaturationBox->setEnabled(true);
}

void LineProfileFitWindow::RangeEdited()
{
    if (nullptr != mLine)
    {
        double *x, *y, *sig;
        int N = mSpektrum->GetLineFitData(x, y, sig, MinFreqEdit->text().toDouble(), MaxFreqEdit->text().toDouble());
        mLine->setData(x, y, sig, N);
        mSpektrum->Changed();
    }
}

void LineProfileFitWindow::RangeChanged(const double newMinE, const double newMaxE)
{
    LineBox->setCurrentIndex(LineBox->count() - 1);
    MinFreqEdit->setText(QString::number(newMinE, 'f', 6));
    MaxFreqEdit->setText(QString::number(newMaxE, 'f', 6));
}
