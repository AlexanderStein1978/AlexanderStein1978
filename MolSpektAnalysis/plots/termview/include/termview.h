//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TERMVIEW_H
#define TERMVIEW_H

#include "tablewidgetwindow.h"

#include <QWidget>

class TermTable;


class TermView : public TableWidgetWindow
{
	Q_OBJECT
public:
    TermView(TermTable *Term);
    virtual ~TermView();

private slots:
	void ShowData();

private:
	double ****Data;
	int vMax, JMax, NComp, NIso;
	TermTable *TT;
};

#endif
