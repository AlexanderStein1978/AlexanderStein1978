//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2026
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
		QPointF position = E->globalPosition();
		emit RightClicked(QPoint(static_cast<int>(position.x()), static_cast<int>(position.y())));
		E->accept();
	}
	else QLineEdit::mousePressEvent(E);
}
