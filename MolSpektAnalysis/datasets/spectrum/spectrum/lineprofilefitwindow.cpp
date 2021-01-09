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

#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QGridLayout>
#include <QIntValidator>
#include <QDoubleValidator>


LineProfileFitWindow::LineProfileFitWindow(MainWindow* mw, Spektrum *spect, Gaussian *line) : LineWindowBase(mw, spect, line),
    MaxIterationEdit(new QLineEdit("100", this)), MinImprovementEdit(new QLineEdit("0.01", this)), MinFreqEdit(new QLineEdit(this)), MaxFreqEdit(new QLineEdit(this)),
    PerformFitButton(new QPushButton("Run fit", this)), ResultSigmaLabel(new QLabel(this))
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
    L->addWidget(PerformFitButton, 5, 0, 1, 2);
    L->addWidget(ResultSigmaLabel, 5, 2, 1, 2);
    MaxIterationEdit->setValidator(new QIntValidator(1, 1000000, MaxIterationEdit));
    MinImprovementEdit->setValidator(new QDoubleValidator(1e-10, 1e10, 0, MinImprovementEdit));
    MinFreqEdit->setValidator(new QDoubleValidator(0.001, 1e10, 0, MinFreqEdit));
    MaxFreqEdit->setValidator(new QDoubleValidator(0.001, 1e10, 0, MaxFreqEdit));
    connect(SpektrumBox, SIGNAL(currentTextChanged(QString)), this, SLOT(SpektrumChanged(QString)));
    connect(LineBox, SIGNAL(currentIndexChanged(int)), this, SLOT(LineChanged(int)));
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
            mSpektrum->FitGaussianLineProfile(lineIndex, MaxIterationEdit->text().toDouble(), MinImprovementEdit->text().toDouble(), MinFreqEdit->text().toDouble(),
                                              MaxFreqEdit->text().toDouble())));
        if (lineIndex >= LineBox->count()) SpektrumChanged(mSpektrum->getFName());
    }
}

void LineProfileFitWindow::connectSpectrum()
{
    if (nullptr != mSpektrum)
    {
        LineBox->addItem("new");
        connect(mSpektrum, SIGNAL(NumberOfFittedLinesChanged()), this, SLOT(NumberOfLinesChanged()));
        connect(mSpektrum, SIGNAL(propertiesChanged()), this, SLOT(UpdateSigma()));
    }
}

void LineProfileFitWindow::disconnectSpectrum()
{
    if (nullptr != mSpektrum)
    {
        disconnect(mSpektrum, SIGNAL(NumberOfFittedLinesChanged()), this, SLOT(NumberOfLinesChanged()));
        disconnect(mSpektrum, SIGNAL(propertiesChanged()), this, SLOT(UpdateSigma()));
    }
}

void LineProfileFitWindow::UpdateSigma()
{
    if (nullptr != mLine) ResultSigmaLabel->setText("Sigma: " + QString::number(mLine->GetSigma()));
}
