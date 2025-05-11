#pragma once

#include "sounddrawwindow.h"


class FrequencyWindow : public SoundDrawWindow
{
    Q_OBJECT
public:
    FrequencyWindow(SoundRecordAndDrawControl *const control, SoundMainWindow *const MW, const int sampleRate);

private slots:
    void BackTransform();
    void CopyAllDataToNewWindow();

private:
    void WriteToFile() override {};
};
