#include "oscillatordataviewer.h"
#include "soundmainwindow.h"
#include "soundwindow.h"
#include "utils.h"
#include "windowselectdialog.h"

#include <QLineEdit>
#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QDoubleValidator>
#include <QComboBox>
#include <QKeyEvent>
#include <QElapsedTimer>

#include <cmath>


OscillatorDataViewer::OscillatorDataViewer(SoundMainWindow* MW, const OscillatorArray::Results& data, const QString& filename, SoundWindow* soundWindow)
    : DiagWindow(SimpleDiagWindow, MW, "Data files (*.dat)", ".dat", 1), mData(data), mLabels(), mSoundWindow(soundWindow), mTimeEdit(new QLineEdit("0.0", this))
    , mStepSizeEdit(new QLineEdit(QString::number(mData.time[1] - mData.time[0], 'g', 5), this)), mLabelBox(new QComboBox(this))
    , mModeButton(new QPushButton("Show", this)), mKeyTime(nullptr), mTimeIndex(0), mHalfDeltaT(0.5 * (data.time[1] - data.time[0]))
{
    setWindowTitle("Oscillator Data Viewer " + filename);
    setUnits("Frequency [Hz]", "Energy [arbitrary units]");
    QPushButton* increaseButton = new QPushButton("+", this), *decreaseButton = new QPushButton("-", this);
    SpektrumLayout->addWidget(new QLabel("Time [s]:", this), 0, 0);
    SpektrumLayout->addWidget(decreaseButton, 0, 1);
    SpektrumLayout->addWidget(mTimeEdit, 0, 2);
    SpektrumLayout->addWidget(increaseButton, 0, 3);
    SpektrumLayout->addWidget(new QLabel("step size:", this), 0, 4);
    SpektrumLayout->addWidget(mStepSizeEdit, 0, 5);
    SpektrumLayout->addWidget(new QLabel("Label:", this), 0, 6);
    SpektrumLayout->addWidget(mLabelBox, 0, 7);
    mModeButton->setCheckable(true);
    SpektrumLayout->addWidget(mModeButton, 0, 8);
    mTimeEdit->setValidator(new QDoubleValidator(data.time[0], data.time[data.numTimeSteps - 1], 10, mTimeEdit));
    mLabelBox->setEditable(false);
    connect(increaseButton, SIGNAL(clicked()), this, SLOT(IncreaseTime()));
    connect(decreaseButton, SIGNAL(clicked()), this, SLOT(DecreaseTime()));
    connect(mTimeEdit, SIGNAL(textChanged(const QString&)), this, SLOT(TimeChanged(const QString&)));
    connect(mLabelBox, SIGNAL(currentIndexChanged(int)), this, SLOT(LabelChanged(int)));
    connect(mModeButton, SIGNAL(clicked(bool)), this, SLOT(ModeChanged(bool)));
    connect(Bild, SIGNAL(KeyPressed(QKeyEvent*)), this, SLOT(KeyPressed(QKeyEvent*)));
}

OscillatorDataViewer::~OscillatorDataViewer()
{
    if (nullptr != mKeyTime) delete mKeyTime;
}

void OscillatorDataViewer::DecreaseTime()
{
    mTimeEdit->setText(QString::number(mTimeEdit->text().toDouble() - mStepSizeEdit->text().toDouble(), 'f', 10));
}

void OscillatorDataViewer::IncreaseTime()
{
    mTimeEdit->setText(QString::number(mTimeEdit->text().toDouble() + mStepSizeEdit->text().toDouble(), 'f', 10));
}

void OscillatorDataViewer::LabelChanged(int index)
{
    mTimeEdit->setText(QString::number(mLabels[index].rect.center().x(), 'f', 10));
}

void OscillatorDataViewer::TimeChanged(const QString& value)
{
    double newTime = value.toDouble();
    while (mTimeIndex < mData.numTimeSteps - 1 && mData.time[mTimeIndex] + mHalfDeltaT < newTime) ++mTimeIndex;
    while (mTimeIndex > 0 && mData.time[mTimeIndex] - mHalfDeltaT > newTime) --mTimeIndex;
    double** currentData = Create(OscillatorArray::NumOscillators, 2);
    for (int n=0; n < OscillatorArray::NumOscillators; ++n)
    {
        currentData[n][0] = mData.frequency[n];
        currentData[n][1] = mData.data[mTimeIndex][n];
    }
    setData(currentData, OscillatorArray::NumOscillators);
}

void OscillatorDataViewer::setLabels(const std::vector<SoundDrawWindow::Label>& labels)
{
    mLabels = labels;
    std::sort(mLabels.begin(), mLabels.end(), [](const SoundDrawWindow::Label& a, const SoundDrawWindow::Label& b)
    {
        return a.rect.center().rx() < b.rect.center().rx();
    });
    mLabelBox->clear();
    for (SoundDrawWindow::Label label : mLabels) mLabelBox->addItem(label.phoneme);
}

void OscillatorDataViewer::KeyPressed(QKeyEvent* K)
{
    QString text = K->text();
    if (text.isEmpty()) return;
    if (mModeButton->isChecked())
    {
        NameSelectionDialog dialog;
        dialog.SetText(text);
        if (dialog.exec() == QDialog::Rejected) return;
        mSoundWindow->addLabel(dialog.GetName(), mData.time[mTimeIndex]);
    }
    else if (mLabelBox->count() > 0)
    {
        if (nullptr == mKeyTime)
        {
            mKeyTime = new QElapsedTimer;
            mKeyTime->start();
            mKeyText = text;
        }
        else if (mKeyTime->restart() < 2000) mKeyText += text;
        else mKeyText = text;
        for (int i = mLabelBox->currentIndex() + 1; i != mLabelBox->currentIndex(); ++i)
        {
            if (i == mLabelBox->count())
            {
                if (mLabelBox->currentIndex() == 0) break;
                i=0;
            }
            if (mLabelBox->itemText(i).compare(mKeyText, Qt::CaseInsensitive) == 0)
            {
                mLabelBox->setCurrentIndex(i);
                break;
            }
        }
    }
}

void OscillatorDataViewer::ModeChanged(bool checked)
{
    if (checked) mModeButton->setText("Set");
    else mModeButton->setText("Show");
}
