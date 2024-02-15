#pragma once

#include "sounddrawwindow.h"


class FrequencyWindow;


class SoundWindow : public SoundDrawWindow
{
    Q_OBJECT

public:
    SoundWindow(SoundRecordAndDrawControl *const control, const QString& filename, const int sampleRate);
    ~SoundWindow();

private slots:
    void Play();
    void FFTActTriggered(bool checked);

private:
    void showFFT() override;
    void WriteToFile() override;
    int getSoundData(float** data);
    int getSoundData(double** data);

    QComboBox* mOutputDeviceBox;
    QAudioOutput* mAudioOutput;
    FrequencyWindow* mFFTWindow = nullptr;
};
