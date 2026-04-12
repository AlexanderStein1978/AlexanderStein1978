//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "Combobox.h"

#include <QMouseEvent>


ComboBox::ComboBox()
{
	Number = 0;
}

ComboBox::ComboBox(QWidget *parent) : QComboBox(parent)
{
	Number = 0;
}

void ComboBox::setNumber(int N)
{
	Number = N;
}

void ComboBox::mouseDoubleClickEvent(QMouseEvent *E)
{
	E->accept();
	emit DoubleClicked(Number);
}

void ComboBox::mousePressEvent(QMouseEvent *E)
{
	if (E->button() == Qt::RightButton)
	{
		QPointF position = E->globalPosition();
		emit RightClicked(QPoint(static_cast<int>(position.x()), static_cast<int>(position.y())));
		E->accept();
	}
	else QComboBox::mousePressEvent(E);
}
