//
// C++ Implementation: TableItem
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "TableItem.h"

#include <QMouseEvent>


TableItem::TableItem() 
{
}

TableItem::TableItem(QString Text) : QLineEdit(Text)
{
}

void TableItem::mousePressEvent(QMouseEvent *E)
{
	if (E->button() == Qt::RightButton)
	{
		emit RightClicked(E->globalPos());
		E->accept();
	}
	else QLineEdit::mousePressEvent(E);
}
