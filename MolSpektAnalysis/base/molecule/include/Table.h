//
// C++ Interface: Table
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_H
#define TABLE_H 


#include <QTableWidget>


class Table : public QTableWidget
{
	Q_OBJECT
			
	public:
		Table(int NR, int NC, QWidget *parent);
	
	protected:
		virtual void keyPressEvent(QKeyEvent *E);
		 
	signals:
		void RowDeleted(bool* delR);
};

#endif
