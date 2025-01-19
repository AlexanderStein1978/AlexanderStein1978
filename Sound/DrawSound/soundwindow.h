#pragma once

#include "sounddrawwindow.h"


class FrequencyWindow;
class QAction;


class SoundWindow : public SoundDrawWindow
{
    Q_OBJECT

    struct IntensityMax
    {
        int F = 0;
        double I = 0.0;
    };

public:
    SoundWindow(SoundRecordAndDrawControl *const control, const QString& filename, const int sampleRate);
    ~SoundWindow();

    void setData(double **Data, int numRows) override;

private slots:
    void play();
    void setFastAssignmentMode(bool enable);
    void continuePlaying();
    void FFTActTriggered(bool checked);
    void AddLabel();
    void LoadLabels();
    void SaveLabels();
    void Delete();
    void WriteAnnInput();
    void ReadAndVerifyAnnOutput();
    void ApplyBoxFilter();
    void ApplyDiffMaxTransfo();
    void mouseLeftClicked(QPoint *Position) override;
    void keyPressed(QKeyEvent *K);

private:
    enum PlayState {PSPlayOnce, PSPlayContinuously, PSStopPlaying};

    void addLabel(const QString name);
    void startPlaying();
    void closeEvent(QCloseEvent *i_event) override;
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    void showFFT() override;
    void WriteToFile() override;
    void CreateAnnInput(double ** data, int &fftLength);
    void ShowPopupMenu(const QPoint& point) override;
    void calcMinLabelWidth();
    void getFFTData(const FFTSelection selection, const int labelIndex, int& N, double **&realData, double **&imaginaryData) const;
    int getSoundData(float** data);
    int getSoundData(double** data, const int labelIndex, FFTSelection fftSelection) const;
    QString predictLabelFilename();
    void analyzeData(double ** const Data, const int numRows);

    QComboBox* mOutputDeviceBox;
    QAudioOutput* mAudioOutput;
    QIODevice* mAudioInputDevice;
    FrequencyWindow* mFFTWindow = nullptr;
    QString mFilename, mLabelFilename, mKeyText;
    const QString mLabelOrderFilename;
    QAction* mAddLabelAct, *mSaveLabelsAct, *mDeleteAct;
    PlayState mPlayState;
    double mMinLabelWidth;
};
