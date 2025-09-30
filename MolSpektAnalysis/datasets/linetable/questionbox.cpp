//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "questionbox.h"

#include <QPushButton>


QuestionBox::QuestionBox(QString Caption, QString Question, QWidget* parent)
			: QMessageBox(QMessageBox::Question, Caption, Question, 
						  QMessageBox::NoButton, parent)
{
	Result = -1;
	NButtons = 0;
}

void QuestionBox::addButton(QString Text)
{
	QMessageBox::ButtonRole R;
	switch (NButtons)
	{
		case 0:
			R = QMessageBox::AcceptRole;
			break;
		case 1:
			R = QMessageBox::NoRole;
			break;
		default:
			R = QMessageBox::RejectRole;
			break;
	}
	QPushButton *B = QMessageBox::addButton(Text, R);
	switch (NButtons)
	{
		case 0:
			connect(B, SIGNAL(clicked()), this, SLOT(b0()));
			break;
		case 1:
			connect(B, SIGNAL(clicked()), this, SLOT(b1()));
			break;
		case 2:
			connect(B, SIGNAL(clicked()), this, SLOT(b2()));
			break;
		case 3:
			connect(B, SIGNAL(clicked()), this, SLOT(b3()));
			break;
		default:
			connect(B, SIGNAL(clicked()), this, SLOT(b4()));
			break;
	}
	NButtons++;
}

void QuestionBox::b0()
{
	Result = 0;
}

void QuestionBox::b1()
{
	Result = 1;
}

void QuestionBox::b2()
{
	Result = 2;
}

void QuestionBox::b3()
{
	Result = 3;
}

void QuestionBox::b4()
{
	Result = 4;
}

int QuestionBox::getResult()
{
	return Result;
}
