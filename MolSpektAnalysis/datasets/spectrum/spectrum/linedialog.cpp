//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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


LineDialog::LineDialog(MainWindow *parent, Spektrum *spect, LineProfile *line)
    : LineWindowBase(parent, spect, line), IntensityEdit(new QLineEdit(this))
    , CenterFreqEdit(new QLineEdit(this)), WidthEdit(new QLineEdit(this)), OffsetEdit(new QLineEdit(this))
    , SubtractButton(new QPushButton(this)), DeleteButton(new QPushButton("Remove line", this))
    , HideSaturationButton(new QPushButton("Hide saturation", this)), DataRangeLabel(new QLabel(this))
{
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
    L->setRowMinimumHeight(5, 20);
    QGridLayout *ButtonLayout = new QGridLayout;
    ButtonLayout->addWidget(SubtractButton, 0, 0);
    ButtonLayout->addWidget(DeleteButton, 0, 1);
    ButtonLayout->addWidget(HideSaturationButton, 0, 2);
    L->addLayout(ButtonLayout, 6, 0, 1, 4);
    connect(SpektrumBox, SIGNAL(currentTextChanged(QString)), this, SLOT(SpektrumChanged(QString)));
    connect(LineBox, SIGNAL(currentIndexChanged(int)), this, SLOT(LineChanged(int)));
    connect(SubtractButton, SIGNAL(clicked()), this, SLOT(SubtractLine()));
    connect(DeleteButton, SIGNAL(clicked()), this, SLOT(DeleteLine()));
    connect(HideSaturationButton, SIGNAL(clicked()), this, SLOT(HideSaturation()));
    connect(IntensityEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    connect(CenterFreqEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    connect(WidthEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    connect(OffsetEdit, SIGNAL(editingFinished()), this, SLOT(UpdateLine()));
    focusInEvent(nullptr);
}

void LineDialog::connectSpectrum()
{
    if (nullptr != mSpektrum)
    {
        connect(mSpektrum, SIGNAL(propertiesChanged()), this, SLOT(lineChanged()));
        connect(mSpektrum, SIGNAL(NumberOfFittedLinesChanged()), this, SLOT(NumberOfLinesChanged()));
    }
}

void LineDialog::disconnectSpectrum()
{
    if (nullptr != mSpektrum)
    {
        disconnect(mSpektrum, SIGNAL(propertiesChanged()), this, SLOT(lineChanged()));
        disconnect(mSpektrum, SIGNAL(NumberOfFittedLinesChanged()), this, SLOT(NumberOfLinesChanged()));
    }
}

void LineDialog::lineChanged()
{
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
        DeleteButton->setEnabled(true);
        SubtractButton->setEnabled(true);
        HideSaturationButton->setEnabled(mLine->isWithSaturation());
    }
    else
    {
        CenterFreqEdit->clear();
        IntensityEdit->clear();
        OffsetEdit->clear();
        WidthEdit->clear();
        DataRangeLabel->clear();
        DeleteButton->setEnabled(false);
        SubtractButton->setEnabled(false);
        HideSaturationButton->setEnabled(false);
    }
    CenterFreqEdit->blockSignals(false);
    IntensityEdit->blockSignals(false);
    OffsetEdit->blockSignals(false);
    WidthEdit->blockSignals(false);
    updateSubtractButton();
}

void LineDialog::UpdateLine()
{
    modifyLine([this]()
    {
        mLine->SetValues(IntensityEdit->text().toDouble(), CenterFreqEdit->text().toDouble(), WidthEdit->text().toDouble(), OffsetEdit->text().toDouble());
    });
}

void LineDialog::updateSubtractButton()
{
    if (nullptr != mLine)
    {
        if (mLine->isLineSubtracted()) SubtractButton->setText("Readd line");
        else SubtractButton->setText("Subtract line");
    }
}

void LineDialog::SubtractLine()
{
    if (nullptr != mSpektrum && nullptr != mLine)
    {
        mSpektrum->SubtractFittedLine(LineBox->currentIndex(), !mLine->isLineSubtracted());
        updateSubtractButton();
    }
}

void LineDialog::DeleteLine()
{
    if (nullptr != mSpektrum && nullptr != mLine) mSpektrum->RemoveFittedLine(LineBox->currentIndex());
}

void LineDialog::selectLine(const int index)
{
    LineBox->setCurrentIndex(index);
}

void LineDialog::HideSaturation()
{
    modifyLine([this]()
    {
        mLine->setHideSaturation(!mLine->isSaturationHidden());
    });
}
