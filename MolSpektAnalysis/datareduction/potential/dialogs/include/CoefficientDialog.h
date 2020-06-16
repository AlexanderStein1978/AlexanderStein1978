//
// C++ Interface: CoefficientDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#ifndef COEFFICIENTDIALOG_H
#define COEFFICIENTDIALOG_H


#include <QDialog>

class QLabel;
class QListWidget;
class QLineEdit;


class CoefficientDialog : public QDialog
{
	Q_OBJECT
	
	public:
		CoefficientDialog(QWidget *parent, QStringList &Coefficients);
		QStringList getResults();
		
	private slots:
		void add();
		void remove();
		
	private:
		QListWidget *CList;
		QLineEdit *Edit;
		QLabel *Label;
		QPushButton *Add, *Remove, *Cancel, *OK;
};

#endif
