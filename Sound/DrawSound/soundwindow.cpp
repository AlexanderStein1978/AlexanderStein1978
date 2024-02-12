#include "soundwindow.h"

#include "utils.h"
#include "fit.h"

#include <QAudioOutput>
#include <QAction>
#include <QMenu>


SoundWindow::SoundWindow(const QString& filename, const int sampleRate) : SoundDrawWindow(sampleRate, 1), mOutputDeviceBox(new QComboBox(this)), mAudioOutput(nullptr)
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioOutput);
    for (QAudioDeviceInfo info : deviceList) mOutputDeviceBox->addItem(info.deviceName());
    QGridLayout *Layout = new QGridLayout;
    Layout->addWidget(new QLabel("Sound output device:", this), 0, 0);
    Layout->addWidget(mOutputDeviceBox, 0, 1);
    SpektrumLayout->addLayout(Layout, 0, 0, 1, 9);
    setUnits("time [s]", "intensity");
    setWindowTitle("Draw sound: " + filename);
    QAction *playAct = new QAction("Play", this), *FFTAct = new QAction("FFT", this);
    FFTAct->setCheckable(true);
    mPopupMenu->addAction(playAct);
    mPopupMenu->addAction(FFTAct);
    connect(playAct, SIGNAL(triggered()), this, SLOT(Play()));
    connect(FFTAct, SIGNAL(toggled(bool)), this, SLOT(FFTActTriggered(bool)));
}


SoundWindow::~SoundWindow() noexcept
{
    if (nullptr != mAudioOutput) delete mAudioOutput;
}

void SoundWindow::Play()
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioOutput);
    QAudioFormat format;
    format.setSampleRate(mSampleRate);
    format.setChannelCount(1);
    format.setSampleSize(32);
    format.setCodec("audio/pcm");
    format.setByteOrder(QAudioFormat::LittleEndian);
    format.setSampleType(QAudioFormat::Float);
    if (nullptr != mAudioOutput) delete mAudioOutput;
    mAudioOutput = new QAudioOutput(deviceList[mOutputDeviceBox->currentIndex()], format, this);
    QIODevice* inputDevice = mAudioOutput->start();
    float* data;
    int length = getSoundData(&data);
    inputDevice->write(reinterpret_cast<char*>(data), 4*length);
}

void SoundWindow::showFFT()
{
    double* data;
    int length = getSoundData(&data), FFTLength = length + 1;
    double **realFFTData = Create(FFTLength, 2), **imaginaryFFTData = Create(FFTLength, 2);
    calcFFT(data, length, 1.0 / mSampleRate, realFFTData, imaginaryFFTData);
    if (nullptr == mFFTWindow)
    {
        mFFTWindow = new DiagWindow;
        mFFTWindow->setUnits("Frequency [Hz]", "Intensity");
    }
    else mFFTWindow->clear();
    mFFTWindow->setData(realFFTData, FFTLength);
    mFFTWindow->addData(imaginaryFFTData, FFTLength);
    if (!mFFTWindow->isVisible()) mFFTWindow->show();
}

void SoundWindow::FFTActTriggered(bool checked)
{
    mIsFFT = checked;
    if (checked && nullptr != mSelectionRect)
    {
        mSelectionRect->setWidth(getFFTWidth(mSelectionRect->width()));
        Paint();
    }
}
