//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "boxfilterdialog.h"

#include <QLineEdit>
#include <QLabel>
#include <QGridLayout>
#include <QPushButton>


BoxFilterDialog::BoxFilterDialog(QWidget* parent) : QDialog(parent), mFilterRadiusEdit(new QLineEdit(this)), mOKButton(new QPushButton("OK", this)), mCancelButton(new QPushButton("Cancel", this))
{
    QGridLayout* L = new QGridLayout(this);
    L->addWidget(new QLabel("Filter radius [samples]:", this), 0, 0);
    L->addWidget(mFilterRadiusEdit, 0, 1);
    L->setRowMinimumHeight(1, 20);
    L->addWidget(mOKButton, 2, 0);
    L->addWidget(mCancelButton, 2, 1);
    connect(mOKButton, SIGNAL(clicked()), this, SLOT(accept()));
    connect(mCancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

int BoxFilterDialog::getResult() const
{
    return mFilterRadiusEdit->text().toInt();
}
