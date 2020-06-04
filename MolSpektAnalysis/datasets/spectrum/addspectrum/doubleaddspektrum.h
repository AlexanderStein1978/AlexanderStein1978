//
// C++ Interface: DoubleAddSpectrum
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2009 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef DOUBLEADDSPEKTRUM_H
#define DOUBLEADDSPEKTRUM_H


#include "Spektrum.h"


class DoubleAddSpectrum : public Spektrum
{
	Q_OBJECT
	
public:
	DoubleAddSpectrum(AddSpectrum *addSpectrum, MainWindow *MW);
    virtual ~DoubleAddSpectrum();
	double getDiff();
	double getSNR();
	
private slots:
	void optimize();
	
	inline void showPMenu(QPoint *P)
	{
		PMenu->popup(*P);
	}
	
private:
	AddSpectrum *addSpect;
	QMenu *PMenu;
};

#endif
