//
// C++ Interface: SourceOffsetDialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SOURCEOFFSETDIALOG_H
#define SOURCEOFFSETDIALOG_H


#include <QDialog> 

class QStringList;
class QLineEdit;

class MainWindow;


class SourceOffsetDialog : public QDialog
{
	Q_OBJECT
	
public:
	SourceOffsetDialog(QStringList &Names, double *Offsets, MainWindow *parent = 0);
	
private slots:
	void ChangeIndex(int Index);
	void Save();
	
private:
	int lastIndex;
	double *offsets;
	QLineEdit *OffsetE;
};

#endif
