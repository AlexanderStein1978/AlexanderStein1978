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



LineDialog::LineDialog(MainWindow *parent, Spektrum *spect, Gaussian *line) : LineWindowBase(parent, spect, line), IntensityEdit(new QLineEdit(this)),
    CenterFreqEdit(new QLineEdit(this)), WidthEdit(new QLineEdit(this)), OffsetEdit(new QLineEdit(this)), DataRangeLabel(new QLabel(this))
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
}



void LineDialog::LineChanged(const int index)
{
    LineWindowBase::LineChanged(index);
    CenterFreqEdit->blockSignals(true);
    IntensityEdit->blockSignals(true);
    OffsetEdit->blockSignals(true);
    WidthEdit->blockSignals(true);
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
