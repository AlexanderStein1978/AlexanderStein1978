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

private:
    QLineEdit* mNameEdit;
    QPushButton *mOkButton, *mCancelButton;
};
