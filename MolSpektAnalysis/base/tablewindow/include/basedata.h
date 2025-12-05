//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#ifndef BASEDATA_H
#define BASEDATA_H


struct BaseData
{
    int progressionNumber = -1, Js = -1, v = -1, J = -1, isotope = -1;
	double energy = 0.0, uncertainty = 0.0, obsMinusCalc = 0.0;
	QString file;
	const QPixmap* IsoIcon = nullptr;
};

#endif
