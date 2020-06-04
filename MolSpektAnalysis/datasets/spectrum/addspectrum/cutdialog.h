//
// C++ Interface: cutdialog
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef CUTDIALOG_H
#define CUTDIALOG_H

/**
	@author Alexander Stein <AlexanderStein@t-online.de>
*/

#include <QWidget>

class MainWindow;
class Spektrum;

class QTableWidget;
class QComboBox;
class QRadioButton;

class CutDialog : public QWidget
{
	Q_OBJECT
			
	public:
    	CutDialog(MainWindow *MW, Spektrum *Spekt);
    	~CutDialog();
		
	private slots:
		void setValue(double x, double y);
		void SpektrumChanged(QString Name);
		void Apply();
		void SpektrumDeleted(int Index);
		
	private:
		MainWindow *MW;
		Spektrum *Spekt;
		QComboBox *SpektrumBox;
		QTableWidget *Tab;
		QRadioButton *in, *out;
};

#endif
