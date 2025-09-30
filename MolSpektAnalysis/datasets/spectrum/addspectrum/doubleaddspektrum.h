//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef DOUBLEADDSPEKTRUM_H
#define DOUBLEADDSPEKTRUM_H


#include "Spektrum.h"

#include <QMenu>

class AddSpectrum;


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
