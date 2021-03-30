//
// C++ Interface: LineDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2020 - 2020
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
    void lineChanged() override;

private:
    void disconnectSpectrum() override;
    void connectSpectrum() override;
    void updateSubtractButton();
    
    QLineEdit *IntensityEdit, *CenterFreqEdit, *WidthEdit, *OffsetEdit;
    QPushButton *SubtractButton, *DeleteButton;
    QLabel *DataRangeLabel;
};

#endif // LINEDIALOG_H
