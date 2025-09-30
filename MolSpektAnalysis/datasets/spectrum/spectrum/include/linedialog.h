//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef LINEDIALOG_H
#define LINEDIALOG_H


#include "linewindowbase.h"

class QLineEdit;
class QLabel;
class QPushButton;

class LineProfile;


class LineDialog : public LineWindowBase
{
    Q_OBJECT
public:
    LineDialog(MainWindow *parent = nullptr, Spektrum* spect = nullptr, LineProfile* line = nullptr);

    void selectLine(const int index);

private slots:
    void UpdateLine();
    void SubtractLine();
    void DeleteLine();
    void HideSaturation();
    void lineChanged() override;

private:
    void disconnectSpectrum() override;
    void connectSpectrum() override;
    void updateSubtractButton();
    void updateHideSaturationButton();
    
    QLineEdit *IntensityEdit, *CenterFreqEdit, *WidthEdit, *OffsetEdit;
    QPushButton *SubtractButton, *DeleteButton, *HideSaturationButton;
    QLabel *DataRangeLabel;
};

#endif // LINEDIALOG_H
