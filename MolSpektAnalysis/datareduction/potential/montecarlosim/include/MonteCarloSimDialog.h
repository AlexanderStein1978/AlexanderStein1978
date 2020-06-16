//
// C++ Interface: MonteCarloSimDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2007 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef MONTECARLOSIMDIALOG_H
#define MONTECARLOSIMDIALOG_H


#include <QDialog>

class MainWindow;

class QLineEdit;


class MonteCarloSimDialog : public QDialog
{
	Q_OBJECT
	
public:
    MonteCarloSimDialog(MainWindow* parent);
	bool exec(QString &Directory, int &NIterations, int &NParFit, double &UncFact);
	
private slots:
	void showFileDialog();
	
private:
	QLineEdit *DirectoryE, *NItE, *NParFitsE, *UncFactE;
};

#endif
