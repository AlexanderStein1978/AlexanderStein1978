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


#include <QWidget>

class QComboBox;
class QLineEdit;
class QLabel;

class Spektrum;
class Gaussian;
class MainWindow;


class LineDialog : public QWidget
{
    Q_OBJECT
public:
    explicit LineDialog(MainWindow *parent = 0);
    void Update();

private slots:
    void SpektrumChanged(const QString& Name);
    void LineChanged(const int index);
    void UpdateLine();

private:
    Spektrum* mSpektrum;
    Gaussian* mLine;
    MainWindow* MW;

    QComboBox *SpektrumBox, *LineBox;
    QLineEdit *IntensityEdit, *CenterFreqEdit, *WidthEdit, *OffsetEdit;
    QLabel *DataRangeLabel;
};

#endif // LINEDIALOG_H
