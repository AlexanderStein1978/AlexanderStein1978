//
// C++ Interface: TableItem
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLEITEM_H
#define TABLEITEM_H 


#include <QLineEdit>


class TableItem : public QLineEdit
{
	Q_OBJECT
	
	public:
		TableItem();
		TableItem(QString);
	
	protected:
		virtual void mousePressEvent(QMouseEvent *event);
	
	signals:
		void RightClicked(QPoint P);
};

#endif
