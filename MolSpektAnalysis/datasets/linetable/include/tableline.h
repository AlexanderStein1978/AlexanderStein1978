//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLELINE_H
#define TABLELINE_H


#include <QString>


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
