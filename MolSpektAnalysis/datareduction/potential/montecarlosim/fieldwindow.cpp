//
// C++ Implementation: FieldWindow
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "fieldwindow.h"


FieldWindow::FieldWindow(int*** iField, int N1, int N2)
{
	QFontMetricsF FM(font());
	NResults = 0;
	Results = 0;
	Step = FM.width("888") + 4;
	Field = iField;
	NxPos = N1;
	NyPos = N2;
	xPos = yPos = 0;
	int height = N2 * Step, width = N1 * Step;
	setMinimumSize(width, height);
	setMaximumSize(width, height);
}

void FieldWindow::paintEvent(QPaintEvent* PE)
{
    //if (Results == 0) return;
	QFontMetricsF FM(font());
	int n, m, c, gOff, yOff, rOff, TOff = (Step + FM.height()) / 2;
	double gSc, ySc, rSc, mg = 1e9, Mg = -1.0, my = 1e9, My = -1.0, mr = 1e9, Mr = -1.0;
	FitResult *FR;
	QPainter P(this);
	QString Text;
	for (n=0; n < NxPos; n++) for (m=0; m < NyPos; m++) if (Field[n][m][0] >= 0)
	{
		FR = Results + Field[n][m][0];
		if (FR->NBad <= 0)
		{
			if (FR->Sigma > Mg) Mg = FR->Sigma;
			if (FR->Sigma < mg) mg = FR->Sigma; 
		}
		else if (FR->NBadPAL == 0)
		{
			if (FR->FQS_Bad > My) My = FR->FQS_Bad;
			if (FR->FQS_Bad < my) my = FR->FQS_Bad;
		}
		else
		{
			if (FR->FQS_Bad > Mr) Mr = FR->FQS_Bad;
			if (FR->FQS_Bad < mr) mr = FR->FQS_Bad;
		}
	}
	if (mg != Mg && Mg >= 0.0)
	{
		gSc = 100.0 / (Mg - mg);
		gOff = int(-mg * gSc);
	}
	else //if (mg == Mg) 
	{
		gSc = 0.0;
		gOff = 0;
	}
	if (my != My && My >= 0.0)
	{
		ySc = 200.0 / (My - my);
		yOff = int(-my * ySc);
	}
	else //if (my == My)
	{
		ySc = 0.0;
		yOff = 0;
	}
	if (mr != Mr && Mr >= 0.0)
	{
		rSc = 255 / (Mr - mr);
		rOff = 255 + int(mr * rSc);
	}
	else //if (mr == Mr)
	{
		rSc = 0.0;
		rOff = 255;
	}
	P.fillRect(0, 0, width(), height(), QColor(255, 255, 255));
	for (n=0; n < NxPos; n++) for (m=0; m < NyPos; m++) if (Field[n][m][0] >= 0)
	{
		FR = Results + Field[n][m][0];
		if (FR->LRelParDev >= 10.0)
		{
			P.fillRect(n * Step, m * Step, Step, Step, QColor(0, 0, 0));
			P.setPen(QColor(255, 255, 255));
		}
		else if (FR->NBad <= 0)
		{
			P.fillRect(n * Step, m * Step, Step, Step, 
					   QColor(0, 255, gOff + int(FR->Sigma * gSc)));
			P.setPen(QColor(0, 0, 0));
		}
		else if (FR->NBadPAL == 0)
		{
			c = yOff + int(FR->FQS_Bad * ySc);
			P.fillRect(n * Step, m * Step, Step, Step, QColor(c, 255 - c, 255 - c));
			P.setPen(QColor(0, 0, 0));
		}
		else
		{
			P.fillRect(n * Step, m * Step, Step, Step, 
					   QColor(rOff - int(FR->FQS_Bad * rSc), 0, 0));
			P.setPen(QColor(255, 255, 255));
		}
		Text = QString::number(Field[n][m][0]);
		P.drawText(n * Step + (Step - FM.width(Text)) / 2, m * Step + TOff, Text);
	}
	if (xPos != 0)
	{
		P.setPen(QColor(255, 0, 0));
		for (n=0; n < MaxParFits; n++) if (xPos[n] >= 0)
		{
			P.fillRect(xPos[n] * Step, yPos[n] * Step, Step, Step, QColor(0, 0, 255));
			Text = QString::number(n);
			P.drawText(xPos[n] * Step + (Step - FM.width(Text)) / 2, 
					   yPos[n] * Step + TOff, Text);
		}
	}
	QWidget::paintEvent(PE);
}

void FieldWindow::update(FitResult* cResults, int nResults, int* Pos1, int* Pos2, int NParFits)
{
	Results = cResults;
	NResults = nResults;
	xPos = Pos1;
	yPos = Pos2;
	MaxParFits = NParFits;
	QWidget::update();
}
