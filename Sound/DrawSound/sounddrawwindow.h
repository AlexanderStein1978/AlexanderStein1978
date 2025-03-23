#pragma once

#include "DiagWindow.h"


class SoundRecordAndDrawControl;
class QAudioOutput;
class SoundMainWindow;


class SoundDrawWindow : public DiagWindow
{
    Q_OBJECT
public:

    struct Label
    {
        QString phoneme;
        QRectF rect;
        int index = -1;
        bool isSelected = false;
    };

    SoundDrawWindow(SoundRecordAndDrawControl *const control, SoundMainWindow *const MW, const int sampleRate, const int o);
    ~SoundDrawWindow();

private slots:
    void mouseMoved(QMouseEvent *e);
    void mousePressed(QMouseEvent *e);
	void mouseReleased(QMouseEvent *e);

protected slots:
    void SelectionChanged(QRect *MarkedArea);
    virtual void WriteToFile() = 0;
    virtual void mouseLeftClicked(QPoint *Position);

protected:
    enum MouseState{MSOutside, MSInside, MSLeft, MSLTCorner, MSTop, MSTRCorner, MSRight, MSRBCorner, MSBottom, MSBLCorner};
    enum FFTSelection{FSForFTT, FSForSelectedFTT, FSNotForFTT, FSDependsOnState};
    enum Mode {MFastLabeling, MNormal};

    void ensureMouseShape(const Qt::CursorShape shape);
    std::pair<int, MouseState> getBestMouseState(const QPoint& point);
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;
    int getFFTLength(const int inputLength) const;
    double getFFTWidth(const double inputWidth) const;
    int getSoundDataRange(int& xStart, int& xStop, const int labelIndex = -1, const FFTSelection fftSelection = FSDependsOnState) const;
    virtual void showFFT() {}
    void ShowPopupMenu(const QPoint& point) override;
    QRect getLabelTextRect(const Label& label);
    int estimateLabelIndex(const QString& phoneme);

    std::vector<Label> mLabels;
    QStringList mLabelOrder;
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
    double mSelectedFFtSize = 0.0;
    Mode mMode;

private:
    static const int D=5;

    void updateLabelRectSelections();
    MouseState isCloseToCorner(const QRectF& rect, const QPoint& point) const;
    MouseState isCloseToWall(const QRectF& rect, const QPoint& point) const;
    bool isInsideRect(const QRectF& rect, const QPoint& point) const;
    bool arePointsClose(const QPointF& pointF, const QPoint& point) const;

    int movingRect;
};
