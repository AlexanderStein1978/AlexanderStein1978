//
// C++ Implementation: LineWindowBase
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "linewindowbase.h"
#include "MainWindow.h"
#include "Spektrum.h"

#include <QComboBox>


LineWindowBase::LineWindowBase(MainWindow *mw, Spektrum *spect, Gaussian *line) : QWidget(mw), SpektrumBox(new QComboBox(this)), LineBox(new QComboBox(this)),
    MW(mw), mSpektrum(spect), mLine(line)
{
    focusInEvent(nullptr);
}

void LineWindowBase::focusInEvent(QFocusEvent *event)
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
    if (nullptr != event) QWidget::focusInEvent(event);
}

void LineWindowBase::SpektrumChanged(const QString &Name)
{
    disconnectSpectrum();
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
                const Gaussian* curLine = nullptr;
                while (curIndex > 0 && curLine != mLine) curLine = mSpektrum->GetFittedLine(--curIndex);
                if (curLine == mLine) LineBox->setCurrentIndex(curIndex);
            }
        }
    }
    LineBox->blockSignals(false);
    connectSpectrum();
    LineChanged(curIndex);
}

void LineWindowBase::LineChanged(const int index)
{
    mLine = (nullptr != mSpektrum && index < mSpektrum->GetNumFittedLines() ? mSpektrum->GetFittedLine(index) : nullptr);
    lineChanged();
}
