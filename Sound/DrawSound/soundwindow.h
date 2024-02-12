#pragma once

#include "sounddrawwindow.h"


class SoundWindow : public SoundDrawWindow
{
    Q_OBJECT

public:
    SoundWindow(const QString& filename, const int sampleRate);
    ~SoundWindow();

private slots:
    void Play();
    void FFTActTriggered(bool checked);

private:
    void showFFT() override;

    QComboBox* mOutputDeviceBox;
    QAudioOutput* mAudioOutput;
    DiagWindow* mFFTWindow = nullptr;
};
