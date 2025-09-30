//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "inputdialog.h"
#include <qlabel.h>
#include <qpushbutton.h>

InputDialog::InputDialog(const double minPeakSize, const double Toleranz, 
			 QWidget *parent, const char */*name*/) : QDialog(parent)
{
    setMinimumSize(200, 120);
    setMaximumSize(200, 120);
    Label1 = new QLabel("Min Peak Size:", this);
    Label1->setGeometry(10, 10, 100, 20);
    Label2 = new QLabel("Suchtoleranz:", this);
    Label2->setGeometry(10, 40, 100, 20);
    minPeak = new QLineEdit(QString::number(minPeakSize), this);
    minPeak->setGeometry(120, 10, 70, 20);
    toleranz = new QLineEdit(QString::number(Toleranz), this);
    toleranz->setGeometry(120, 40, 70, 20);
    QPushButton *OK = new QPushButton("OK", this);
    OK->setGeometry(10, 90, 80, 20);
    QPushButton *Cancel = new QPushButton("Cancel", this);
    Cancel->setGeometry(110, 90, 80, 20);
    connect(OK, SIGNAL(clicked()), this, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), this, SLOT(reject()));
}

InputDialog::~InputDialog()
{

}

void InputDialog::setLabel(QString L1, QString L2)
{
	Label1->setText(L1);
	Label2->setText(L2);
}

void InputDialog::getData(double &V1, double &V2)
{
	V1 = minPeak->text().toDouble();
	V2 = toleranz->text().toDouble();
}

void InputDialog::getData(int &v1, int &v2)
{
	v1 = minPeak->text().toInt();
	v2 = toleranz->text().toInt();
}

double InputDialog::GetMinPeak() const
{
    return minPeak->text().toDouble();
}

double InputDialog::GetToleranz() const
{
    return toleranz->text().toDouble();
}
