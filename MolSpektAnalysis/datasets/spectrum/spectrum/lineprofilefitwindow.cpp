//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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


LineProfileFitWindow::LineProfileFitWindow(MainWindow* mw, Spektrum *spect, LineProfile *line) : LineWindowBase(mw, spect, line),
    MaxIterationEdit(new QLineEdit("100", this)), MinImprovementEdit(new QLineEdit("0.01", this)), MinFreqEdit(new QLineEdit(this)),
    MaxFreqEdit(new QLineEdit(this)), LineTypeBox(new QComboBox(this)), PerformFitButton(new QPushButton("Run fit", this)),
    ResultSigmaLabel(new QLabel(this)), mLineDialog(nullptr), mWithSaturationBox(new QCheckBox("Consider saturation", this))
{
    QGridLayout *L = new QGridLayout(this);
    L->addWidget(new QLabel("Spectrum:", this), 0, 0);
    L->addWidget(SpektrumBox, 0, 1);
    L->addWidget(new QLabel("Line:", this), 0, 2);
    L->addWidget(LineBox, 0, 3);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Profile type:", this), 2, 0);
    L->addWidget(LineTypeBox, 2, 1);
    for (int t = LineProfile::GaussianType; t <= LineProfile::LorentzianType; ++t)
        LineTypeBox->addItem(LineProfile::getProfileTypeName(static_cast<LineProfile::LineProfileType>(t)));
    LineTypeBox->setEditable(false);
    L->addWidget(mWithSaturationBox, 2, 2, 1, 2);
    L->addWidget(new QLabel("Max iterations:", this), 3, 0);
    L->addWidget(MaxIterationEdit, 3, 1);
    L->addWidget(new QLabel("Min improvement:", this), 3, 2);
    L->addWidget(MinImprovementEdit, 3, 3);
    L->addWidget(new QLabel("Min frequency:", this), 4, 0);
    L->addWidget(MinFreqEdit, 4, 1);
    L->addWidget(new QLabel("Max frequency:", this), 4, 2);
    L->addWidget(MaxFreqEdit, 4, 3);
    L->setRowMinimumHeight(5, 20);
    L->addWidget(PerformFitButton, 6, 0, 1, 2);
    L->addWidget(ResultSigmaLabel, 6, 2, 1, 2);
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
        modifyLine([&lineIndex,this]()
        {
            ResultSigmaLabel->setText("Sigma: " + QString::number(
                mSpektrum->FitLineProfile(lineIndex, static_cast<LineProfile::LineProfileType>(LineTypeBox->currentIndex()),
                                          mWithSaturationBox->isChecked(), MaxIterationEdit->text().toDouble(),
                                          MinImprovementEdit->text().toDouble(), MinFreqEdit->text().toDouble(),
                                          MaxFreqEdit->text().toDouble())));
        });
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
    }
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
