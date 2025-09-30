//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

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
