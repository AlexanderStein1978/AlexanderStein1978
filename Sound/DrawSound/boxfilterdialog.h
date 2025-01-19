#pragma once

#include <QDialog>


class QLineEdit;


class BoxFilterDialog : public QDialog
{
    Q_OBJECT
public:
    BoxFilterDialog(QWidget *parent = nullptr);

    int getResult() const;

private:
    QLineEdit* mFilterRadiusEdit;
    QPushButton* mOKButton, *mCancelButton;
};
