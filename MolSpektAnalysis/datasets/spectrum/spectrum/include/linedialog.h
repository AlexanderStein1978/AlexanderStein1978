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


class LineDialog : public LineWindowBase
{
    Q_OBJECT
public:
    explicit LineDialog(MainWindow *parent = nullptr, Spektrum* spect = nullptr, Gaussian* line = nullptr);

private slots:
    void UpdateLine();
    void lineChanged() override;

private:
    void disconnectSpectrum() override;
    void connectSpectrum() override;
    
    QLineEdit *IntensityEdit, *CenterFreqEdit, *WidthEdit, *OffsetEdit;
    QLabel *DataRangeLabel;
};

#endif // LINEDIALOG_H
