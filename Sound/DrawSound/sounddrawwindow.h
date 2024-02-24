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
    void mouseLeftClicked(QPoint *Position);

protected slots:
    virtual void WriteToFile() = 0;

protected:
    enum MouseState{MSOutside, MSInside, MSLeft, MSLTCorner, MSTop, MSTRCorner, MSRight, MSRBCorner, MSBottom, MSBLCorner};

    struct Label
    {
        QString phoneme;
        QRectF rect;
        bool isSelected = false;
    };

    void ensureMouseShape(const Qt::CursorShape shape);
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    int getFFTLength(const int inputLength);
    double getFFTWidth(const double inputWidth);
    int getSoundDataRange(int& xStart, int& xStop);
    virtual void showFFT() {}
    void ShowPopupMenu(const QPoint& point) override;
    QRect getLabelTextRect(const Label& label);

    std::vector<Label> mLabels;
    QMenu *mPopupMenu;
    QRectF *mSelectionRect = nullptr;
    int mSampleRate;
    bool mIsFFT = false;
    SoundRecordAndDrawControl *const mControl;
    QRectF mMoveStartRect;
    MouseState mMouseState = MSOutside, mMoveState = MSOutside;
    QPoint mMoveMouseStartPoint;
    QFont mLabelFont;
    QRectF* mMovingRect = nullptr;

private:
    static const int D=5;

    void updateLabelRectSelections();
    std::pair<int, MouseState> getBestMouseState(const QPoint& point);
    MouseState isCloseToCorner(const QRectF& rect, const QPoint& point) const;
    MouseState isCloseToWall(const QRectF& rect, const QPoint& point) const;
    bool isInsideRect(const QRectF& rect, const QPoint& point) const;
    bool arePointsClose(const QPointF& pointF, const QPoint& point) const;


    int movingRect;
};
