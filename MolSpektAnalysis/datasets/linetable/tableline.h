//
// C++ Interface: tableline
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2017
//
// Copyright: See README file that comes with this source code
//
//


#ifndef TABLELINE_H
#define TABLELINE_H


class LineTable;
class ElState;


struct TableLine
{
	int vs, Js, vss, Jss, Iso, PN, Row, SourceN, FC, nDig;
	double WN, err, dev, DevR, SNR, WeightFact;
	QString File, SourceName;
	LineTable *LTab;
    ElState *State;
	bool isTE, isViewn, isSelected;
};

#endif
