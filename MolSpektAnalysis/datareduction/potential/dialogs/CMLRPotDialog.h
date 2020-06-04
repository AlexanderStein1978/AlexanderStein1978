//
// C++ Interface: CMLRPotDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef CMLRPOTDIALOG_H
#define CMLRPOTDIALOG_H


#include <QDialog>


class CMLRPotDialog : public QDialog
{
	Q_OBJECT
	
public:
    CMLRPotDialog(QStringList *CVL, QWidget* parent = 0);
	
	QLineEdit *pE, *qE, *RrefE, *mRE, *MRE, *NPE, *NCE;
	QListWidget *LRC;
	
private slots:
	void addC();
	void delC();
	void setV();
	
private:
	QComboBox *LRCB;
	QPushButton *addB, *delB, *setB, *OK, *Cancel;
	QStringList *CVL;
	QLineEdit *CValE;
};

#endif
