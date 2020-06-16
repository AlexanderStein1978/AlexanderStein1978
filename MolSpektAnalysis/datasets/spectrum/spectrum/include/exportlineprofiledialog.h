//
// C++ Interface: ExportLineProfileDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef EXPORTLINEPROFILEDIALOG_H
#define EXPORTLINEPROFILEDIALOG_H


#include <QWidget>

class Spektrum;
class MainWindow;

class QLineEdit;
class QComboBox;
class QPushButton;


class ExportLineProfileDialog : public QWidget
{
    Q_OBJECT
public:
    ExportLineProfileDialog(Spektrum* const i_spectrum, MainWindow* i_MW);

private slots:
    void LineChanged(int i_index);
    void Write();
    void SelectFileName();

signals:
    void closeThis();

private:
    Spektrum *m_spectrum;
    QLineEdit *m_FileNameEdit, *m_eStartEdit, *m_eEndEdit, *m_eStepEdit;
    QComboBox *m_lineSelectBox;
    QPushButton *m_OKButton, *m_CancelButton, *m_selectFileButton;

};

#endif
