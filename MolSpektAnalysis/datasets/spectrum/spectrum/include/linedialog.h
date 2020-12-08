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
    virtual void LineChanged(const int index) override;

private:
    void UpdateLine();

    QLineEdit *IntensityEdit, *CenterFreqEdit, *WidthEdit, *OffsetEdit;
    QLabel *DataRangeLabel;
};

#endif // LINEDIALOG_H
