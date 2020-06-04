//
// C++ Interface: SpectList
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#ifndef SPEKTLIST_H
#define SPEKTLIST_H


#include "tablewindow.h"


class SpectList : public TableWindow
{
public:
	SpectList(MainWindow *MW = 0);
	void AutoSLP();
};

#endif
