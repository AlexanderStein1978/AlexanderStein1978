//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
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
