//
// C++ Implementation: LineDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "linedialog.h"
#include "Spektrum.h"
#include "gaussian.h"
#include "MainWindow.h"

#include <QComboBox>
#include <QLineEdit>
#include <QLabel>
#include <QGridLayout>
#include <QDoubleValidator>



LineDialog::LineDialog(MainWindow *parent) : QWidget(parent), mSpektrum(nullptr), mLine(nullptr), MW(parent), SpektrumBox(new QComboBox(this)),
    LineBox(new QComboBox(this)), IntensityEdit(new QLineEdit(this)), CenterFreqEdit(new QLineEdit(this)), WidthEdit(new QLineEdit(this)),
    OffsetEdit(new QLineEdit(this)), DataRangeLabel(new QLabel(this))
{
    SpektrumBox->setEditable(false);
    LineBox->setEditable(false);
    IntensityEdit->setValidator(new QDoubleValidator(IntensityEdit));
    CenterFreqEdit->setValidator(new QDoubleValidator(CenterFreqEdit));
    WidthEdit->setValidator(new QDoubleValidator(WidthEdit));
    OffsetEdit->setValidator(new QDoubleValidator(OffsetEdit));
    QGridLayout *L(new QGridLayout(this));
    L->addWidget(new QLabel("Spectrum: ", this), 0, 0);
    L->addWidget(SpektrumBox, 0, 1);
    L->addWidget(new QLabel("Line:", this), 0, 2);
    L->addWidget(LineBox, 0, 3);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Intensity:", this), 2, 0);
    L->addWidget(IntensityEdit, 2, 1);
    L->addWidget(new QLabel("Center frequency:", this), 2, 2);
    L->addWidget(CenterFreqEdit, 2, 3);
    L->addWidget(new QLabel("Width:", this), 3, 0);
    L->addWidget(WidthEdit, 3, 1);
    L->addWidget(new QLabel("Offset:", this), 3, 2);
    L->addWidget(OffsetEdit, 3, 3);
    L->addWidget(DataRangeLabel, 4, 0, 1, 4);
    connect(SpektrumBox, SIGNAL(currentTextChanged(QString)), this, SLOT(SpektrumChanged(QString)));
    connect(LineBox, SIGNAL(currentIndexChanged(int)), this, SLOT(LineChanged(int)));
    connect(IntensityEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    connect(CenterFreqEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    connect(WidthEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    connect(OffsetEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    Update();
}

void LineDialog::Update()
{
    SpektrumBox->blockSignals(true);
    SpektrumBox->clear();
    if (nullptr != MW)
    {
        int NSpectra = MW->getNumSpectra();
        for (int n=0; n < NSpectra; ++n)
        {
            Spektrum* curSpect = MW->getSpectrum(n);
            if (curSpect->GetNumFittedLines() > 0)
            {
                SpektrumBox->addItem(curSpect->getFName());
                if (mSpektrum == curSpect) SpektrumBox->setCurrentIndex(SpektrumBox->count() - 1);
            }
        }
    }
    SpektrumBox->blockSignals(false);
    SpektrumChanged(SpektrumBox->currentText());
}

void LineDialog::SpektrumChanged(const QString &Name)
{
    LineBox->blockSignals(true);
    int curIndex = LineBox->currentIndex();
    LineBox->clear();
    if (nullptr != MW)
    {
        int NSpectra = MW->getNumSpectra();
        for (int n=0; n < NSpectra; ++n)
        {
            Spektrum* curSpect = MW->getSpectrum(n);
            if (curSpect->getFName() == Name)
            {
                mSpektrum = curSpect;
                break;
            }
        }
        if (nullptr != mSpektrum)
        {
            int nLines = mSpektrum->GetNumFittedLines();
            for (int n=0; n < nLines; ++n) LineBox->addItem(QString::number(n));
            if (nullptr != mLine)
            {
                if (curIndex > nLines) curIndex = nLines;
                else ++curIndex;
                Gaussian* curLine = nullptr;
                while (curIndex > 0 && curLine != mLine) curLine = mSpektrum->GetFittedLine(--curIndex);
                if (curLine == mLine) LineBox->setCurrentIndex(curIndex);
            }
        }
    }
    LineBox->blockSignals(false);
    LineChanged(curIndex);
}

void LineDialog::LineChanged(const int index)
{
    CenterFreqEdit->blockSignals(true);
    IntensityEdit->blockSignals(true);
    OffsetEdit->blockSignals(true);
    WidthEdit->blockSignals(true);
    mLine = (nullptr != mSpektrum && index < mSpektrum->GetNumFittedLines() ? mSpektrum->GetFittedLine(index) : nullptr);
    if (nullptr != mLine)
    {
        double CenterFreq, Intensity, Offset, Width, FStart, FEnd;
        mLine->GetValues(Intensity, CenterFreq, Width, Offset);
        mLine->GetDataRange(FStart, FEnd);
        CenterFreqEdit->setText(QString::number(CenterFreq, 'f'));
        IntensityEdit->setText(QString::number(Intensity, 'f'));
        OffsetEdit->setText(QString::number(Offset, 'f'));
        WidthEdit->setText(QString::number(Width, 'f'));
        DataRangeLabel->setText(QString("Data range: %1 cm-1 to %2 cm-1").arg(FStart, 0, 'f', 6).arg(FEnd, 0, 'f', 6));
    }
    else
    {
        CenterFreqEdit->clear();
        IntensityEdit->clear();
        OffsetEdit->clear();
        WidthEdit->clear();
        DataRangeLabel->clear();
    }
    CenterFreqEdit->blockSignals(false);
    IntensityEdit->blockSignals(false);
    OffsetEdit->blockSignals(false);
    WidthEdit->blockSignals(false);
}

void LineDialog::UpdateLine()
{
    if (nullptr != mLine)
    {
        mLine->SetValues(IntensityEdit->text().toDouble(), CenterFreqEdit->text().toDouble(), WidthEdit->text().toDouble(), OffsetEdit->text().toDouble());
        mSpektrum->Changed();
    }
}
