#pragma once

#include "DiagWindow.h"


class QAudioOutput;


class SoundDrawWindow : public DiagWindow
{
    Q_OBJECT
public:
    SoundDrawWindow(const QString& filename, const int sampleRate);
    ~SoundDrawWindow();

private slots:
    void SelectionChanged(QRect *MarkedArea);
    void mouseMoved(QMouseEvent *e);
    void mousePressed(QMouseEvent *e);
	void mouseReleased(QMouseEvent *e);
    void Play();
    void WriteToFile();

private:
    enum MouseState{MSOutside, MSInside, MSLeft, MSLTCorner, MSTop, MSTRCorner, MSRight, MSRBCorner, MSBottom, MSBLCorner};

    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    void ShowPopupMenu(const QPoint& point) override;
    void ensureMouseShape(const Qt::CursorShape shape);
    int getSoundData(float** data);

    QRectF *mSelectionRect = nullptr, mMoveStartRect;
    MouseState mMouseState = MSOutside, mMoveState = MSOutside;
    QPoint mMoveMouseStartPoint;
    QMenu *mPopupMenu;
    QComboBox* mOutputDeviceBox;
    QAudioOutput* mAudioOutput;
    int mSampleRate;
};
