#pragma once

#include "sounddrawwindow.h"


class FrequencyWindow;
class QAction;


class SoundWindow : public SoundDrawWindow
{
    Q_OBJECT

public:
    SoundWindow(SoundRecordAndDrawControl *const control, const QString& filename, const int sampleRate);
    ~SoundWindow();

private slots:
    void Play();
    void FFTActTriggered(bool checked);
    void AddLabel();
    void LoadLabels();
    void SaveLabels();

private:


    void closeEvent(QCloseEvent *i_event) override;
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    void showFFT() override;
    void WriteToFile() override;
    void ShowPopupMenu(const QPoint& point) override;
    int getSoundData(float** data);
    int getSoundData(double** data);
    QString predictLabelFilename();

    QComboBox* mOutputDeviceBox;
    QAudioOutput* mAudioOutput;
    FrequencyWindow* mFFTWindow = nullptr;
    QString mFilename;
    QString mLabelFilename;
    QAction* mAddLabelAct, *mSaveLabelsAct;
};
