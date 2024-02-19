#pragma once

#include "DiagWindow.h"


class SoundRecordAndDrawControl;
class QAudioOutput;


class SoundDrawWindow : public DiagWindow
{
    Q_OBJECT
public:
    SoundDrawWindow(SoundRecordAndDrawControl *const control, const int sampleRate, const int o);
    ~SoundDrawWindow();

private slots:
    void SelectionChanged(QRect *MarkedArea);
    void mouseMoved(QMouseEvent *e);
    void mousePressed(QMouseEvent *e);
	void mouseReleased(QMouseEvent *e);
    virtual void WriteToFile() = 0;

protected:
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    int getFFTLength(const int inputLength);
    double getFFTWidth(const double inputWidth);
    int getSoundDataRange(int& xStart, int& xStop);
    virtual void showFFT() {}
    void ShowPopupMenu(const QPoint& point) override;

    QMenu *mPopupMenu;
    QRectF *mSelectionRect = nullptr;
    int mSampleRate;
    bool mIsFFT = false;
    SoundRecordAndDrawControl *const mControl;

private:
    enum MouseState{MSOutside, MSInside, MSLeft, MSLTCorner, MSTop, MSTRCorner, MSRight, MSRBCorner, MSBottom, MSBLCorner};

    void ensureMouseShape(const Qt::CursorShape shape);

    QRectF mMoveStartRect;
    MouseState mMouseState = MSOutside, mMoveState = MSOutside;
    QPoint mMoveMouseStartPoint;
};
