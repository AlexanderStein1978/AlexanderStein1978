#include "soundwindow.h"
#include "frequencywindow.h"
#include "datensatz.h"
#include "utils.h"
#include "fit.h"

#include <QAudioOutput>
#include <QAction>
#include <QMenu>
#include <QFileDialog>


SoundWindow::SoundWindow(SoundRecordAndDrawControl *const control, const QString& filename, const int sampleRate) : SoundDrawWindow(control, sampleRate, 1), mOutputDeviceBox(new QComboBox(this)),
    mAudioOutput(nullptr)
{
    QList<QAudioDeviceInfo> deviceList = QAudioDeviceInfo::availableDevices(QAudio::AudioOutput);
    for (QAudioDeviceInfo info : deviceList) mOutputDeviceBox->addItem(info.deviceName());
    mOutputDeviceBox->setEditable(false);
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


SoundWindow::~SoundWindow()
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
    delete[] data;
}

int SoundWindow::getSoundData(float ** data)
{
    int xStart, xStop, n, i, rLength = getSoundDataRange(xStart, xStop);
    *data = new float[rLength];
    for (n = xStart, i=0; n <= xStop; ++n, ++i) (*data)[i] = static_cast<float>(Daten->GetValue(n, 1));
    return rLength;
}

int SoundWindow::getSoundData(double ** data)
{
    int xStart, xStop, n, i, rLength = getSoundDataRange(xStart, xStop);
    *data = new double[rLength];
    for (n = xStart, i=0; n <= xStop; ++n, ++i) (*data)[i] = Daten->GetValue(n, 1);
    return rLength;
}

void SoundWindow::showFFT()
{
    double* data;
    int length = getSoundData(&data), FFTLength = length + 1;
    double **realFFTData = Create(FFTLength, 2), **imaginaryFFTData = Create(FFTLength, 2);
    calcFFT(data, length, 1.0 / mSampleRate, realFFTData, imaginaryFFTData);
    if (nullptr == mFFTWindow) mFFTWindow = new FrequencyWindow(mControl, mSampleRate);
    else mFFTWindow->clear();
    mFFTWindow->setData(realFFTData, FFTLength);
    mFFTWindow->addData(imaginaryFFTData, FFTLength);
    if (!mFFTWindow->isVisible()) mFFTWindow->show();
    Destroy(realFFTData, FFTLength);
    Destroy(imaginaryFFTData, FFTLength);
    delete[] data;
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

void SoundWindow::WriteToFile()
{
    QString filename = QFileDialog::getSaveFileName(this, "Select filename to save", DATA_DIRECTORY  "/Chars");
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    float* data;
    int length = getSoundData(&data);
    file.write(reinterpret_cast<char*>(data), 4*length);
    delete[] data;
}
