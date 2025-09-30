//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TERMENERGY_H
#define TERMENERGY_H


class ElState;


struct TermEnergy
{
	int v, J, Iso, PN, FC, nDig, Row;
	double E, err, dev, DevR, WeightFactor;
	bool ef;
	QString File;
    ElState *State;
};

#endif
