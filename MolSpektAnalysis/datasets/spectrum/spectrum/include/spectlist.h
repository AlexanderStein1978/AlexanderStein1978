//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef SPEKTLIST_H
#define SPEKTLIST_H


#include "tablewidgetwindow.h"


class SpectList : public TableWidgetWindow
{
	Q_OBJECT
public:
	SpectList(MainWindow *MW = 0);
	void AutoSLP();
};

#endif
