//
// C++ Implementation: ComboBox
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "ComboBox.h"

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
		emit RightClicked(E->globalPos());
		E->accept();
	}
	else QComboBox::mousePressEvent(E);
}