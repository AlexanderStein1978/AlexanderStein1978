#pragma once

#include "sounddrawwindow.h"


class FrequencyWindow : public SoundDrawWindow
{
    Q_OBJECT
public:
    FrequencyWindow(SoundRecordAndDrawControl *const control, const int sampleRate);

private slots:
    void BackTransform();
    void CopyAllDataToNewWindow();

private:
    void WriteToFile() override {};
};
