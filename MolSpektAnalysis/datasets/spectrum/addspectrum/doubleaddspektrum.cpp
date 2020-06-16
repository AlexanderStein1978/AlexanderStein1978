//
// C++ Implementation: DoubleAddSpectrum
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "doubleaddspektrum.h"
#include "datensatz.h"
#include "addspektrum.h"


DoubleAddSpectrum::DoubleAddSpectrum(AddSpectrum* addSpectrum, MainWindow* MW): Spektrum(MW)
{
	addSpect = addSpectrum;
	PMenu = new QMenu;
	QAction *OptimizeAct = new QAction("&Optimize", this);
	PMenu->addAction(OptimizeAct); 
	connect(Bild, SIGNAL(RightClicked(QPoint*)), this, SLOT(showPMenu(QPoint*)));
	connect(OptimizeAct, SIGNAL(triggered()), this, SLOT(optimize()));
}

DoubleAddSpectrum::~DoubleAddSpectrum()
{
	delete PMenu;
}

double DoubleAddSpectrum::getDiff()
{
	int n = Daten->GetDSL() / 2;
	double y1, y2, y3, x, y;
	while (Daten->GetValue(n, 0) > 0) n--;
	while (Daten->GetValue(n, 0) < 0) n++;
	if (Daten->GetValue(n-1, 1) > Daten->GetValue(n+1, 1))
	{
		for (y1 = Daten->GetValue(n, 1), y2 = Daten->GetValue(--n, 1);
				 y2 < (y3 = Daten->GetValue(--n, 1)); y1 = y2, y2 = y3) ;
		ParabInterpol(Daten->GetValue(n+2, 0), y1, Daten->GetValue(n+1, 0), y2,
					  Daten->GetValue(n, 0), y3, x, y);
	}
	else
	{
		for (y1 = Daten->GetValue(n, 1), y2 = Daten->GetValue(++n, 1);
				 y2 < (y3 = Daten->GetValue(++n, 1)); y1 = y2, y2 = y3) ;
		ParabInterpol(Daten->GetValue(n-2, 0), y1, Daten->GetValue(n-1, 0), y2,
					  Daten->GetValue(n, 0), y3, x, y);
	}
	return x;
}

void DoubleAddSpectrum::optimize()
{
	addSpect->maximizeDAdd(getDiff());
}

double DoubleAddSpectrum::getSNR()
{
	editFind();
	Marker *Line = getMarker(0.0);
	return (Line != 0 ? Line->SNR : 0.0);
}
