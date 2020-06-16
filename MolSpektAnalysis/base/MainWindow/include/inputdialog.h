//
// C++ Interface: InputDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef INPUTDIALOG_H
#define INPUTDIALOG_H


#include <qdialog.h>
#include <qlineedit.h>

class QLabel;

class InputDialog : public QDialog
{
public:
    InputDialog(const double minPeakSize, const double Toleranz, 
		QWidget *parent = NULL, const char *name = NULL);
    ~InputDialog();
    void setLabel(QString Label1, QString Label2);
	void getData(double &V1, double &V2);
	void getData(int &v1, int &v2);
	double GetMinPeak() const;
    double GetToleranz() const;
private:
    QLineEdit *minPeak, *toleranz;
	QLabel *Label1, *Label2;
};

#endif
