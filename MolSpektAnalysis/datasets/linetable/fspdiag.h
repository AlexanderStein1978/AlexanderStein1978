//
// C++ Interface: FSPDiag
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef FSPDIAG_H
#define FSPDIAG_H


#include <QWidget>

class QComboBox;
class QPushButton;

class MainWindow;


class FSPDiag : public QWidget
{
	Q_OBJECT
	
	public:
		FSPDiag(MainWindow *MW, Progression P);
		~FSPDiag();
		
	private slots:
		void search();
		
	private:
		QComboBox *LTBox;
		QPushButton *SearchB, *CloseB;
		Progression *P;
		MainWindow *MW;
};

#endif
