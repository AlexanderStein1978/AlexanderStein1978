//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once

#include <QDialog>


class DiagWindow;
class QComboBox;
class QLineEdit;


class WindowSelectDialog : public QDialog
{
    Q_OBJECT

public:
    WindowSelectDialog(const std::vector<DiagWindow*>& windows, QWidget *parent = nullptr);

    int GetSelection();

private:
    QComboBox* mSelectionBox;
    QPushButton *mOkButton, *mCancelButton;
};


class NameSelectionDialog : public QDialog
{
    Q_OBJECT

public:
    NameSelectionDialog(QWidget *parent = nullptr);

    QString GetName();
    void SetText(const QString& text);

private:
    QLineEdit* mNameEdit;
    QPushButton *mOkButton, *mCancelButton;
};
